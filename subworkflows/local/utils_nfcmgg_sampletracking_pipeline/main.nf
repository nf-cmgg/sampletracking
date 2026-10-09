def haplotypeMapToBed(inputFilePath: Path, outputFilePath: Path) {
    outputFilePath.text = ""
    inputFilePath.eachLine { line ->
        // Skip metadata headers (@) and column definition lines (#)
        if (line.startsWith('@') || line.startsWith('#')) {
            return null
        }
        def parts = line.split('\t')
        if (parts.size() >= 3) {
            def chrom: String = parts[0]
            def position: long = parts[1].toLong()
            def name: String = parts[2]
            // BED format: 0-based start, 1-based end
            def chromStart: long = position - 1
            def chromEnd: long = position
            // Write standard 4-column BED (chrom, start, end, name)
            outputFilePath.append("${chrom}\t${chromStart}\t${chromEnd}\t${name}\n")
        }
    }
}

def depthFilter(sample_depth_tsv: Path, snp_depth_tsv: Path, min_covered_sites: int) {
    // read both files, split by tabs
    def sampleDepthMap = sample_depth_tsv.readLines().collectEntries { line ->
        def parts = line.split('\t')
        [(parts[0]+":"+parts[1]): parts[2].toInteger()]
    }
    def snpDepthMap = snp_depth_tsv.readLines().collectEntries { line ->
        def parts = line.split('\t')
        [(parts[0]+":"+parts[1]): parts[2].toInteger()]
    }
    def coveredSites = sampleDepthMap.keySet().intersect(snpDepthMap.keySet())
    println("Covered sites: ${coveredSites}")
    // count the number of sites covered in both files in the same positions
    if (coveredSites.size() >= min_covered_sites) {
        return true
    }
    return false
}
