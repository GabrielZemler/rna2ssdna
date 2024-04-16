# rna2ssdna
rna2ssdna project

The aim of this project is to help to generate ssDNA hairpins and aptamers for molecular dynamics simulation by changing structures obtained from RNAcomposerer etc. (pdb format) and translating them to DNA. All will be tested with the AMBER99sb-ILDN force field in GROMCACS 2022.3.

In the src folder you will find a precompiled executable (rna2ssdna) that you can use. If the executable is not working on your system you can always compile it yourself. Since the code contains the floating-point versions of some mathematical functions like cosf, srqtf etc. in math.h you may have to add "-lm" to use the correct header file. 

You can add the path to the executable to your .bashrc, -zshrc or profile etc. to make the program executable from your working directories. If you use it for the first time, try "rna2ssdna -h" which gives you some helpful info about the possible parameters. A generic usage would be "rna2ssdna -f yourRNA.pdb. This will generate an "output.pdb" and a "rna2ssdna.log" file. The first one contains your new DNA structure, while the second one contains a full log of the calculations, insertions and replacements that were done. 

It is still a work in progress and I will implement more features as time goes on. Please always check with your viewer of choice to see if the changes were made correctly. If you encounter problems, feel free to email me at Gabriel.Zemler@biologie.uni-regensburg.de. If you use this for academic interest, please mention the use of this program in your references.

Good luck!