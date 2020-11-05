# Program Assignment of Course: Combinatorial Methods in Bioinformatics
# Author: 生信1702 邵燕涵 2017317220205
# Software: PyCharm 2020.1.1 (Community Edition)
# Environment: python3.7


一、From the report.pdf you can check all of the result report coming from testing data and download data from http://eggnogdb.embl.de/#/app/downloads. Besides, there are also contains some other analysis and explanations. Take a look at it first and it will save you a lot of time to testing all of the algorithms!


二、From the codes folder you can find two python codes file:
1、linear_gap_penalty.py:
	This file is the algorithm to handle linear gap penalty, including two functions for input sequence and parameter file reading, one function for back-tracking, three functions for DP banded-DP and X-drop and last the one is main function. To testing this algorithm all you need to do is to fill in the path of your input file, parameter file and three output file in the main function area. After that, just run this file so you will get the result file and the program's result report! For more details please check the codes' annotations.
	Notice: Because the input file: input1.txt and input3.txt, parameter file: parameter1.txt and parameter3.txt are nucleotide sequences, and their score for initiating a gap are 0, so I using them as the testing data in main function. But if you want, you can also using the protein sequences input2.txt and parameter2.txt as a testing data, and before that, you should remember to change the file path and change the nucleotide sequence reading function to protein input sequence reading function, because there is a little difference between them.

2、affine_gap_penalty.py:
	This file is the algorithm to handle affine gap penalty, including same parts of linear gap penalty.py. The only difference is this algorithm need to calculate three table V, E and F to handle the affine gap penalty model. To testing this algorithm, all you need to do is to fill in the path of your input file, parameter file and three output file in the main function area. After that, just run this file so you will get the result file and the program's result report! For more details please check the codes' annotations.
	Notice: Because the input file: input2.txt and parameter file: parameter2.txt are protein sequences, and the score for initiating a gap are not 0, so I using them as the testing data in main function. But if you want, you can also using the nucleotide  sequences input1.txt, input3.txt, parameter1.txt and parameter3.txt as a testing data, but remember their score for initiating a gap are 0, maybe you need to change it first, and before that, you should also remember to change the file path and change the protein sequence reading function to nucleotide  input sequence reading function, because there is a little difference between them.


三、From the input folder you can find the testing input and parameter file: 
Input1.txt, Input2.txt, Input3.txt 
parameter1.txt, parameter2.txt, parameter3.txt


四、From the output folder you can find the output file calculated from the testing input file:
out1_DP.txt, out1_band.txt, out1_drop.txt
out2_DP.txt, out2_band.txt, out2_drop.txt
out3_DP.txt, out3_band.txt, out3_drop.txt


五、I also chose some homology protein sequences from http://eggnogdb.embl.de/#/app/downloads for affine gap penalty algorithm's analysis, you can read my report.pdf directly or testing it by yourself. You can find the data under folder analytical_input:
analysis_input1.txt, analysis_input2.txt, analysis_input3.txt, analysis_input4.txt
parameter.txt







