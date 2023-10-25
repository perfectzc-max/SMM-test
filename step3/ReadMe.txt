#1. Identification of UMI families. Command line parameters (positional, whithout keys)
main.py umi <S file> <Chrom (use '-' for None)>

A report is generated for each UMI family:
UMI, Start, End, strang count +, strang count -.

#2. Variant calling. Command line parameters (positional, whithout keys)
main.py vc <S file> <Chrom (use '-' for None)> <G file> <Reference> <SNP database> <G VCF records file>

Result formats are presented in the project description.

The additional parameters are listed in the file Params.ini. Their meaning is clear enough from the names and values.
Some explanations are given below. 

The family includes all reads with a common UMI, the areas [nStart, nEnd) of which intersect at least at maxPosDist positions.
A related pair read1, read2 is united if they contents at least one common position.

Qualified families are determined from the size of each strand +,- (each >= minStrandCount)
or the total size of the family ( >= minFamilyCount, if this parameter is positive).

File names and chromosomes can also be specified by parameters in Params.ini (see examples there). But these parameters can
"overlapped" by command line parameters. You can set up to 6 of them, see the main.py file and comments there.

#3.usage of each py file
3.1 main.py
main.py 用于读取输入的运行参数，包括运行的任务如umi、vc，提供的文件地址。并调用umiExplorer里的方法。
输入对应的参数并调用umiExplorer.umiAnalyze()
3.1.1 UMI
命令：
main.py umi bam/PBMM02_S.gatk.bam -
输入一个样本的bam文件，和对应的CHROM，然后就会直接调用umiExplorer.umiAnalyze()，最后会给出一个UMIfamily的报告，但感觉不是很有必要，可以运行看看，但不是必须的。

3.1.1 VC
命令：
main.py vc <S file> <Chrom (use '-' for None)> <G file> <Reference> <SNP database> <G VCF records file>
这个会输出一个记录了mutation的VCF文件，路径就在bam文件同路径。
3.2 params.ini
记录了各种参数，在运行的时候可以在里面改好文件位置，以便后续运行。如果定义好了文件位置，就不用反复输入同样的参考序列地址。
3.3 params.py
文件定义了各种所需的参数，包括文件、命令、算法处理的具体参数。初始的参数文件params.ini里的参数值将会被调用。
3.4 umiexplorer.py
创建UMIExplorer这个个object。
调用了family.py和utils.py文件里的各种方法。
3.5 family.py
调用了utils里的方法。
创建Read、Family两个object，并定义一系列方法。
3.6 utils.py
创建Utils这个个object。
定义了一些基本方法，详见文件内注释。





