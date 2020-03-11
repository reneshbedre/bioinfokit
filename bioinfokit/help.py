class format:
    def __init__(self):
        pass

    @staticmethod
    def fq_qual_var():
        text = """
            <b>FASTQ quality format detection</b>

            `bioinfokit.analys.format.fq_qual_var(file)`
                
             Parameters:
             ------------
             file : FASTQ file to detect quality format [deafult: None]
            
            Returns:
            Quality format encoding name for FASTQ file
            (Supports only Sanger, Illumina 1.8+ and Illumina  1.3/1.4)

            <a href="https://reneshbedre.github.io/blog/fqqualfmt.html" target="_blank">Working Example</a>
            """
        print(text)