import java.util.zip.GZIPInputStream

/**
 * Contig-compatibility checks between the reference genome (.fai) and the packaged resources a run will use.
 *
 * Pure Groovy: Nextflow loads it from lib/ and main.nf calls ReferenceCheck.run(...) at the top level, so the check
 * runs on the head node before any process is defined.  Only headers are read for VCF, .dict and .interval_list
 * files; BED-like files are streamed once.
 *
 * Rules
 *   dict-style resources (VCF ##contig header, .dict, .interval_list @SQ header):
 *       must share at least one contig with the reference, and every shared contig must have the same length.
 *       Extra contigs on either side are allowed (alt / decoy / HLA differences between GRCh38 builds).
 *   span-style resources (BED, TSV, interval bodies, VCF records when the header has no ##contig lines):
 *       every contig used must exist in the reference (unless allowUnknown) and no coordinate may exceed its length.
 *   --ref_dict: must list exactly the contigs of --ref_fai with the same lengths.
 */
class ReferenceCheck {

    static BufferedReader open(String path) {
        InputStream stream = new FileInputStream(new File(path))
        if (path.endsWith('.gz') || path.endsWith('.bgz')) stream = new GZIPInputStream(stream, 1 << 16)
        return new BufferedReader(new InputStreamReader(stream, 'UTF-8'), 1 << 16)
    }

    /** .fai: name <TAB> length ... */
    static Map<String, Long> fromFai(String path) {
        def m = [:]
        open(path).withReader { r ->
            r.eachLine { line ->
                if (!line.trim()) return
                def f = line.split('\t')
                m[f[0]] = f[1] as Long
            }
        }
        return m
    }

    /** SAM-style header (.dict, .interval_list): @SQ SN:name LN:length */
    static Map<String, Long> fromSqHeader(String path) {
        def m = [:]
        open(path).withReader { r ->
            String line
            while ((line = r.readLine()) != null) {
                if (!line.startsWith('@')) break
                if (!line.startsWith('@SQ')) continue
                String sn = null; Long ln = null
                line.split('\t').each { tag ->
                    if (tag.startsWith('SN:')) sn = tag.substring(3)
                    if (tag.startsWith('LN:')) ln = tag.substring(3) as Long
                }
                if (sn != null && ln != null) m[sn] = ln
            }
        }
        return m
    }

    /** VCF header: ##contig=<ID=name,length=N,...>  ->  name -> length (null if the header has no contig lines) */
    static Map<String, Long> fromVcfHeader(String path) {
        def m = [:]
        boolean any = false
        open(path).withReader { r ->
            String line
            while ((line = r.readLine()) != null) {
                if (!line.startsWith('#')) break
                if (!line.startsWith('##contig=<')) continue
                any = true
                def body = line.substring('##contig=<'.length(), line.lastIndexOf('>'))
                String id = null; Long len = null
                body.split(',').each { kv ->
                    int i = kv.indexOf('=')
                    if (i < 0) return
                    def k = kv.substring(0, i), v = kv.substring(i + 1)
                    if (k == 'ID') id = v
                    if (k == 'length' && v.isLong()) len = v as Long
                }
                if (id != null) m[id] = len
            }
        }
        return any ? m : null
    }

    /** Tabular body (BED, TSV, interval_list body, VCF records): contig -> highest coordinate seen */
    static Map<String, Long> spanFromTable(String path, long maxLines) {
        def m = [:]
        long n = 0
        open(path).withReader { r ->
            String line
            while ((line = r.readLine()) != null) {
                if (!line || line.startsWith('#') || line.startsWith('@') || line.startsWith('track') || line.startsWith('browser')) continue
                def f = line.split('\t')
                if (f.length < 2) continue
                Long end = null
                if (f.length >= 3 && f[2].isLong()) end = f[2] as Long       // BED / interval_list end
                else if (f[1].isLong())              end = f[1] as Long       // 2-column loci / VCF POS
                if (end == null) continue
                def cur = m[f[0]]
                if (cur == null || end > cur) m[f[0]] = end
                if (++n >= maxLines) break
            }
        }
        return m
    }

    static List<String> compareDict(Map<String, Long> ref, Map<String, Long> res, String label, boolean sameSet = false) {
        def problems = []
        if (!res) { problems << "${label}: no contigs found in its header"; return problems }
        def common = res.keySet().intersect(ref.keySet())
        if (!common) {
            problems << "${label}: no contig in common with the reference (resource: ${res.keySet().take(3).join(', ')} ...; reference: ${ref.keySet().take(3).join(', ')} ...)"
            return problems
        }
        def bad = common.findAll { res[it] != null && res[it] != ref[it] }
        if (bad) problems << "${label}: ${bad.size()} shared contig(s) with a different length, e.g. ${bad.take(3).collect { "${it} ${res[it]} vs ${ref[it]} in reference" }.join('; ')}"
        if (sameSet) {
            def onlyRes = res.keySet() - ref.keySet()
            def onlyRef = ref.keySet() - res.keySet()
            if (onlyRes) problems << "${label}: ${onlyRes.size()} contig(s) not in --ref_fai, e.g. ${onlyRes.take(3).join(', ')}"
            if (onlyRef) problems << "${label}: ${onlyRef.size()} --ref_fai contig(s) missing, e.g. ${onlyRef.take(3).join(', ')}"
        }
        return problems
    }

    static List<String> compareSpan(Map<String, Long> ref, Map<String, Long> span, String label, boolean allowUnknown) {
        def problems = []
        if (!span) { problems << "${label}: no records found"; return problems }
        def unknown = span.keySet() - ref.keySet()
        if (unknown.size() == span.size()) {
            problems << "${label}: none of its contigs exist in the reference (e.g. ${unknown.take(3).join(', ')})"
            return problems
        }
        if (unknown && !allowUnknown) problems << "${label}: ${unknown.size()} contig(s) not in the reference, e.g. ${unknown.take(3).join(', ')}"
        def over = span.findAll { c, end -> ref.containsKey(c) && end > ref[c] }
        if (over) problems << "${label}: coordinates beyond the contig length, e.g. ${over.take(3).collect { c, e -> "${c}:${e} > ${ref[c]}" }.join('; ')}"
        return problems
    }

    /**
     * resources: list of maps
     *   [label:, path:, kind: 'dict'|'vcf'|'interval_list'|'bed', severity: 'error'|'warn', allowUnknown: bool, maxLines: long]
     *
     * severity says what the consuming tool does on a contig mismatch:
     *   'error' (default) - the tool aborts (GATK validates sequence dictionaries and -L intervals against the reference)
     *   'warn'            - the tool runs but silently ignores the unmatched regions or mis-annotates (SAGE / mosdepth /
     *                       Strelka / Manta BED lookups, MuSE dbSNP annotation, Delly exclude list)
     * A missing or unreadable file is always an error because every tool opens its resource file.
     *
     * Returns [errors: [...], warnings: [...]]; both empty means compatible.  Empty paths are skipped.
     */
    static Map run(String refFai, List<Map> resources) {
        def errors = [], warnings = []
        Map<String, Long> ref
        try { ref = fromFai(refFai) } catch (Exception e) { return [errors: ["cannot read --ref_fai ${refFai}: ${e.message}".toString()], warnings: []] }
        if (!ref) return [errors: ["--ref_fai ${refFai} lists no contigs".toString()], warnings: []]

        resources.each { res ->
            // empty name (e.g. RN7 has no PON / dbSNP) leaves main.nf with "" or "<database dir>/": nothing to check
            if (!res.path || res.path.toString().endsWith('/')) return
            String path = res.path.toString()
            String label = "${res.label} (${path})"
            def file = new File(path)
            if (file.isDirectory()) return
            if (!file.exists()) { errors << "${label}: file not found".toString(); return }
            long maxLines = (res.maxLines ?: Long.MAX_VALUE) as long
            boolean allowUnknown = res.allowUnknown ?: false
            def found = []
            try {
                switch (res.kind) {
                    case 'dict':
                        found.addAll(compareDict(ref, fromSqHeader(path), label, true)); break
                    case 'interval_list':
                        found.addAll(compareDict(ref, fromSqHeader(path), label))
                        found.addAll(compareSpan(ref, spanFromTable(path, maxLines), label, allowUnknown)); break
                    case 'vcf':
                        def hdr = fromVcfHeader(path)
                        if (hdr != null) found.addAll(compareDict(ref, hdr, label))
                        else found.addAll(compareSpan(ref, spanFromTable(path, Math.min(maxLines, 100000L)), label, allowUnknown))
                        break
                    case 'bed':
                        found.addAll(compareSpan(ref, spanFromTable(path, maxLines), label, allowUnknown)); break
                    default:
                        errors << "${label}: unknown resource kind '${res.kind}'".toString(); return
                }
            } catch (Exception e) {
                errors << "${label}: cannot be read (${e.message})".toString(); return
            }
            if (found) (res.severity == 'warn' ? warnings : errors).addAll(found)
        }
        return [errors: errors.collect { it.toString() }, warnings: warnings.collect { it.toString() }]
    }
}
