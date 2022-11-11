# DualMindex 81.20.0101

Release Date: 13th September, 2022

Author

	~Changyong Yu (Northeastern University in CHINA)

1.Introduction

	DualMindex is a tool，it can construct the index of de bruijn graph (cDBG) of one genome or multiple genomes. The built index can represent the genomes in the data analysis such as large-scale read mapping.

	The k-mers of the genome sequences are the vertex of the graph, the length of the k-mer is an important parameter of the graph. DualMindex limits the length of k-mer less than or equal to 128. In our experience, 128 is long enough for various kinds of aims of data analysis when using cDBG. But also, we can modify the codes to support longer k-mer vertex.

	The DualMindex is composed of two parts of information, 1) the partial FM-index of a minimal sequence g' which is a shortest sequence generating the same DBG with the analyzed genomes; 2) the position list of the minimal edges. In our opinior, for large-scale, especially huge-scale genome or multiple genomes, we recommend FM-index combining with position lists of minimal edges for indexing the position of the edges for reducing the cost of memory. 

	DualMindex is a powerful tool. You can use DualMindex to generate the efficient index of large-scale genome sequences. The main advantage of DualMindex is that it is smaller than most of the available index meanwhile it is fast for index DBG and calculate the position list of any path. In our experiments it can deal with up to 7 human-size genomes with a pc of 160 GB memory on acceptable time cost, the generated index size is no more than 36.5GB.

	There is a shortcoming of DualMindex. DualMindex is designed base on the assumption that the alphbet is {A,C,G,T}. Therefore, DualMindex can just recognize the four letters. It cannot deal with the other characters with recognizing them. If it encounter with other characters, it will treat them as letter 'A'. This surely will result in false positive results. The further verification on the genome sequence will correct and delete the false positive results. It is worth attention for the users of cDBG for this shortcoming. We think it is not vitally important compared with the key problem of how to filter out true negitive results in data analysis such as read mapping. Therefore, it is accecptable.

	Overall, DualMindex can construct an effective index of DBG for multiple large-scale genomes for read mapping. We believe that the tool will improve the performance of data analysis with the bottleneck of large-scale data.

2.Test Data

	DualMindex constructs the index of the de bruijn graph of genome sequences. The genome sequences are .fa format, and we tested DualMindex with the genome sequences used in paper "TwoPaCo: An efficient algorithm to build the compressed de Bruijn graph from many complete genomes". we download the data using website: https://github.com/medvedevgroup/TwoPaCo.

3.Building Notes

Our directory contains only two contents: 1.codes to build DualMindex of de bruijn graph; 2.codes to build FM-index of MNR-path (tool named gMap_fm).

THe directories named #src# and #src_fm# consisted of the code necessary to build the DualMindex and gMap_fm.
To build the DualMindex and gMap_fm, change the directory to src (src_fm) and type

	$ make

After that, you can type commands to run tests using StLiter. StLiter provides 5 methods for building  compressed de bruijn graph of large-scale genomes.

4.Usage Notes

1)O method. Generate the optimal MNR-path of the DBG.

	$ ./gMap -M <method> -kL <kmer_len> -R <ref_path>
	
	O method calcultes effectively the MNR-paths. It output one file, take kL=16 for example, it will output one files, "aa16". It is the file saving the MNR-path with k=16. The NR-paths in the files are connected by char '#'. 

Running test:

	$ ./gMap -M O -kL 16 -R dataset1

	The above commond takes "dataset1" file and "nomLeu2.fa" file as input, it will output file named "aa16". 

2)o method. Generate the optimal MNR-path of the DBG from an MNR-path of the genomes corresponding to longer k-mer.

	$ ./gMap -M <method> -kL <kmer_len> -R <ref_path>
	
	o method calcultes effectively the MNR-paths. Different with 'O' method, it take an MNR-path the genomes with longer k-mer as input. It output one file, take kL=16 for example, it will output one files, "aa16". It is the file saving the MNR-path with k=16. The NR-paths in the files are connected by char '#'. 

Running test:

	$ ./gMap -M O -kL 16 -R aa22

	It will output file named "aa16". The reference sequences are corresponding to the ones that generate the input file "aa22".

3)V method. Generate the bit vector for labelling the minimal edges.

	$ ./gMap -M <method> -kL <kmer_len> -t <thread_num> -R <ref_path> -Vg <bit gap> 

	V method calculates the bit vector for labelling the minimal edges. It outputs two files, one is the bit vector and the other is its index. The parameter "bit gap" must be divisible by 8. "bit gap" is a parameter of generated file Nb**. When the generated file Nb** is used in the following, the parameter used must be the same with the one in this commond. 

Running test:

	$ ./gMap -M V -kL 16 -t 4 -R aa16 -Vg 256
		
	Input files:
		B, C, OCCA and aa16

	Output files:
		vb16 and Nb16

	The above commonds calculate bit vector for labelling minimal edges in DBG with k-mer length 16. 

4)L method. Generate the position lists of the minimal edges.

	 $ ./gMap -M <method> -kL <kmer_len> -Vg <bit gap> -Vb <bv file> -Nb <Nb file> -t <thread num> -R <ref_path>

	L method generate the position list and its index.

Running test:

	$ ./gMap -M L -kL 16 -R dataset7 -Vg 256 -Vb vb16 -Nb Nb16 -t 1

	Input files:
		B,C OCCA,dataset7(and these genome sequences),vb16,Nb16

	Output files:
		Pi16 and Ps16

	The above commonds generate the position lists of minimal edges and it index file. 

5)F method. Locate any path using the generated index.

	$ ./gMap -M <method> -kL <kmer_len> -Vb <bv file> -Nb <Nb file> -Pi <Pi file> -Ps <Ps file> -Vg <bit gap> -R <ref_path>

	F method calcultes the position list of a sub-string of reference. 

Running test:

	$ ./gMap -M F -kL 16 -Vb vb16 -Nb Nb16 -Pi Pi16 -Ps Ps16 -Vg 256 -R <ref_path>

	Input files:
		B,C,OCCA,vb16,Nb16,Pi16,Ps16

	Output files:
		a position list of a sub-string of references

	The above commonds calculate the position list of the input path.

6)B method. Build the FM-index of a string (the alphabet is {A,C,G,T,$})

	$ ./gMap_fm -M <method> -h <string_path> -s <SA_gap> -o <OCC_gap> -l <level> -t <thread_num>

	B method calcultes the compressed FM-index of a string on the alphabet {A,C,G,T,$}. 

Running test:

	$ ./gMap_fm -M B -h aa10 -s 10 -o 50 -l 2 -t 0

	Input files:
		aa10

	Output files:
		B,C,OCCA,SA(not used in DualMindex)

	The above commonds calculate the compressed FM-index of the input string. note that the thread_num is on log number 4^thread_num. the level parameter decide the partition of the suffices, and each partition is calculated in a single loop. This method can efficient reduce the memory cost.

7)T method. test the FM-index of a string (the alphabet is {A,C,G,T,$})

	$ ./gMap_fm -M <method> -h <string_path>

	T method test the compressed FM-index of a string on the alphabet {A,C,G,T,$}. 

Running test:

	$ ./gMap_fm -M T -h aa10

	Input files:
		aa10,B,C,OCCA,SA

	Output files:
		none

	The above commonds test the compressed FM-index of the input string.

8)I method. build a compressed FM-index based on a non-compressed FM-index.

	$ ./gMap_fm -M <method> -s <SA_gap> -o <OCC_gap>

	T method build a new compressed FM-index. 

Running test:

	$ ./gMap_fm -M I -s 10 -o 50

	Input files:
		B,C,OCCA,SA

	Output files:
		B,cC,cOCCA,cSA (with compressing parameter SA_gap=10, OCC_gap=50)

	The above commonds build a new compressed FM-index of the input string.


5.Parameter Settings

The format of a parameter of DualMindex in the command line is a pair of strings, here we denote the pair as (-p, [q]) or (-p,<q>). String p is the name of the parameter. String q is the value of the parameter input in the command line. [q] represents that the parameter is a optional parameter. <q> represents that the parameter is a necessary parameter.

@parameter (-M,<method>)

	Parameter 'M' assigns the specific method, such as 'O','o','V','L','F' and 'C'.

@parameter (-R,<ref_path>)

	Parameter 'r' gives the path of a text file which saves the filenames of the genome files. For example 'ref_path'="dataset1", than dataset1 is a text file which saves the filename "nomLeu2.fa".

@parameter (-kL,<kmer_len>)

	Parameter 'kL' set the length of the k-mer

@parameter (-t,<thread_num>)

	The parameter 't' sets the number of threads. 

@parameter (-Vb,<vb file>)

	Parameter 'Vb' assigns the file path of the bit vector file vb**.

@parameter (-Nb,<Nb file>)

	Parameter 'Nb' assigns the file path of the bit vector file Nb**.

@parameter (-Pi,<Pi file>)

	Parameter 'Pi' assigns the file path of the bit vector file Pi**.

@parameter (-Ps,<Ps file>)

	Parameter 'Ps' assigns the file path of the bit vector file Ps**.

6.Format of cDBG

	-aaxx
		~binary files of MNR-path
		~xx is the k-mer length
		~each Byte is a character of MNR-path.

	-B,C,OCCA
		~binary files of partial FM-index

	-vbxx
		~binary file of bit vector of minimal edges
		~xxx is the k-mer length
		~the i-th bit is 1 if the edge with id i is a minimal edge.
		
	-Nbxx
		~binary file of vbxx index
		~xxx is the k-mer length
		~store uint32_t integers, Nb[i] for labelling the number of '1's in the first i*Vg_gap bits in vbxx. 

	-Psxx
		~binary file of position lists
		~xx is the k-mer length
		~store uint32_t integers for position lists of all minimal edges, sorted by the id of minimal edge.
	-Pixx
		~binary file of index of position lists
		~xx is the k-mer length
 		~store uint32_t integers, Pi[i] for labelling the starting position of postion list of the i-th minimal edge.

7.License

	See LICENSE.txt

8.Contacts

	Please e-mail your feedback at cyyneu@126.com.



