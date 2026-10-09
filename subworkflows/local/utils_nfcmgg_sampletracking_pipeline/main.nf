def haplotypeMapToBed(inputFilePath: Path, outputFilePath: String) {
    def outputFile = new File(outputFilePath)

    outputFile.withWriter('utf-8') { writer ->
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
                writer.writeLine("${chrom}\t${chromStart}\t${chromEnd}\t${name}")
            }
        }
    }
}
