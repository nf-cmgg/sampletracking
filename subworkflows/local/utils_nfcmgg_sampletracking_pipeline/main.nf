def haplotypeMapToBed(String inputFilePath, String outputFilePath) {
    def inputFile = new File(inputFilePath)
    def outputFile = new File(outputFilePath)

    outputFile.withWriter('utf-8') { writer ->
        inputFile.eachLine { line ->
            // Skip metadata headers (@) and column definition lines (#)
            if (line.startsWith('@') || line.startsWith('#')) {
                return
            }

            def parts = line.split('\t')
            if (parts.size() >= 3) {
                String chrom = parts[0]
                long position = parts[1].toLong()
                String name = parts[2]

                // BED format: 0-based start, 1-based end
                long chromStart = position - 1
                long chromEnd = position

                // Write standard 4-column BED (chrom, start, end, name)
                writer.writeLine("${chrom}\t${chromStart}\t${chromEnd}\t${name}")
            }
        }
    }
}
