import os, sys, glob


#EX) ARF4  1       43718   43919   1:43818 999     .       96.0    24.59   22.70   100
TFBS1=['1_43718_43919']
Not_TFBS1=['1_43500_43718']
Not_TFBS2=['1_20000_20500']
TFBS2=['1_55316_55517']
Not_TFBS3=['1_55000_55300']
TFBS3=['1_175739_175940']
TFBS4=['1_337794_337995']
TFBS5=['1_440655_440856']
TFBS6=['1_464133_464334']
Not_TFBS4=['1_464335_465000']


#What I imagine to make input file
#chr_Pos Seq-1 Seq(e.g.A,T,C,G) Seq+1  TFBS //DistanceFromGene//-->later  Depth-3 Depth-2 Depth -1 Depth0 Depth+1 Depth+2 Depth+3 

Outfile = open("ExampleSequence_Matrix","w")
##Maize 


def MakeFastaDic():
	Dic = {}
	infile = open("Chr1_00","r")
	#>1 dna:chromosome chromosome:AGPv4:1:1:307041717:1 REF
	for sLine in infile:
		if sLine.startswith(">"):
			if " " in sLine:
				Chr = sLine.strip().split(" ")[0].replace(">","")
			else:
				Chr = sLine.strip().replace(">","")
			Dic.setdefault(Chr,"")

		else:
			Dic[Chr]+= sLine.strip()
	return Dic
			
def ParseSequence(List,TFBS_Switch):
	for sPos in List:
		sPos_list = sPos.split("_")
		sChr=sPos_list[0]
		nStart=int(sPos_list[1])
		nEnd=int(sPos_list[2])
	
		while (nStart<nEnd):
			sChrPos = sChr+"_"+str(nStart)
			
			Seq_m1=FASTADic[sChr][nStart-1]
			Seq=FASTADic[sChr][nStart]
			Seq_p1=FASTADic[sChr][nStart+1]

			Depth_m3=DepthDic[sChr][nStart-3]
			Depth_m2=DepthDic[sChr][nStart-2]
			Depth_m1=DepthDic[sChr][nStart-1]
			Depth=DepthDic[sChr][nStart]
			Depth_p1=DepthDic[sChr][nStart+1]
			Depth_p2=DepthDic[sChr][nStart+2]
			Depth_p3=DepthDic[sChr][nStart+3]
			
			Outfile.write(sChrPos+"\t"+Seq_m1+"\t"+Seq+"\t"+Seq_p1+"\t")
			Outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(Depth_m3,Depth_m2,Depth_m1,Depth,Depth_p1,Depth_p2,Depth_p3,TFBS_Switch))

			nStart+=1
	
	
def DepthFile():
	Dic = {}
	infile = open("ARF4_Depth_00","r")
	for sLine in infile:
		sList =sLine.strip().split("\t")
		#print(sList)
		Dic.setdefault(sList[0],[]) ## List !! 
		Dic[sList[0]].append(int(sList[2]))
	return Dic
	

FASTADic = MakeFastaDic()
print("Make fasta dic done!")
print(FASTADic.keys())
DepthDic = DepthFile()
print("Make Depth dic done!")
Outfile.write("chr_Pos\tSeq_m1\tSeq\tSeq_p\tDepth_m3\tDepth_m2\tDepth_m1\tDepth\tDepth_p1\tDepth_p2\tDepth_p3\tTFBS\n")
ParseSequence(TFBS1,"TFBS")
ParseSequence(TFBS2,"TFBS")
ParseSequence(TFBS3,"TFBS")
ParseSequence(TFBS4,"TFBS")
ParseSequence(TFBS5,"TFBS")
ParseSequence(TFBS6,"TFBS")

ParseSequence(Not_TFBS1,"Not_TFBS")
ParseSequence(Not_TFBS2,"Not_TFBS")
ParseSequence(Not_TFBS3,"Not_TFBS")
ParseSequence(Not_TFBS4,"Not_TFBS")


Outfile.close()
