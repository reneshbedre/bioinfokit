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

class stat:
    def __init__(self):
        pass

    @staticmethod
    def lin_reg():
        text = """
            <b>Linear regression analysis</b>

            `bioinfokit.analys.stat.linearreg(file)`

             Parameters:
             ------------
             df: Pandas dataframe object
             x : Name of column having independent X variables [list][default:None]
             y : Name of column having dependent Y variables [list][default:None]

            Returns:
            Regression analysis summary

            <a href="https://reneshbedre.github.io/blog/linearreg.html" target="_blank">Working Example</a>
            """
        print(text)

    @staticmethod
    def regplot():
        text = """
                <b>Regression plot</b>

                `bioinfokit.visuz.stat.regplot(df, x, y, yhat, dim, colordot, colorline, r, ar, dotsize, markerdot, linewidth, 
                    valphaline, valphadot)`

                 Parameters:
                 ------------
                 df        : Pandas dataframe object
                 x         : Name of column having independent X variables [string][default:None]
                 y         : Name of column having dependent Y variables [string][default:None]
                 yhat      : Name of column having predicted response of Y variable (y_hat) from regression [string][default:None]
                 dim       : Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
                 r         : Figure resolution in dpi [int][default: 300]
                 ar        : Rotation of X-axis labels [float][default: 0]
                 dotsize   : The size of the dots in the plot [float][default: 6]
                 markerdot : Shape of the dot marker. See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
                valphaline : Transparency of regression line on plot [float (between 0 and 1)][default: 1]
                valphadot  : Transparency of dots on plot [float (between 0 and 1)][default: 1]
                linewidth  : Width of regression line [float][default: 1]

                Returns:

                Regression plot image in same directory (reg_plot.png)

                <a href="https://reneshbedre.github.io/blog/linearreg.html" target="_blank">Working Example</a>
                """
        print(text)