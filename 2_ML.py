import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import datasets
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn import tree
from sklearn.model_selection import cross_val_score
from sklearn import preprocessing
from sklearn.feature_extraction.text import CountVectorizer
#module load Python/3.7.0-foss-2018a

def Convert(Input):
	if Input == "A":
   	    	Output = 1
	elif Input == "T":
		Output = 2
	elif Input == "G":
               	Output = 3
	elif Input == "C":
               	Output = 4
	else:
		print(Input)
	return Output

def SingleInput(InfileName):
	X=[]
	y=[]
	infile = open(InfileName,"r")
	infile.readline()
	for sLine in infile:
		sList = sLine.strip().split("\t")
		#subList = [Convert(sList[1]),Convert(sList[2]),Convert(sList[3])]
		#subList = [int(sList[4]),int(sList[5]),int(sList[6]),int(sList[7]),int(sList[8]),int(sList[9])]
		subList = [sList[1],sList[2],sList[3]]
		X.append(subList)
		y.append(sList[11])

	#one-hot encoding
	
	enc = preprocessing.OneHotEncoder()
	enc.fit(X)	
	NewX =  enc.transform(X).toarray()
	
	return NewX,y	

def SingleInput_DAP_Depth(InfileName):
	X=[]
	y=[]
	infile = open(InfileName,"r")
	infile.readline()
	for sLine in infile:
		sList = sLine.strip().split("\t")
		subList = [int(sList[4]),int(sList[5]),int(sList[6]),int(sList[7]),int(sList[8]),int(sList[9])]
		X.append(subList)
		y.append(sList[11])
	return X,y



def MultipleInput(InfileName):
	X=[]
	y=[]
	infile = open(InfileName,"r")
	infile.readline()
	for sLine in infile:
		sList = sLine.strip().split("\t")
		Seq_windo=sList[1]
		Seq_windo_new =""
		#for sSeq in Seq_windo:
			#Seq_windo_new+=str(Convert(sSeq))
		#print(Seq_windo_new)
		#X.append(int(Seq_windo_new))
		X.append([Seq_windo])
		                vectorizer = CountVectorizer()
                print(subList)
                Temp = vectorizer.fit_transform(subList)
                print(Temp.toarray())

		y.append(sList[2])
	infile.close()
	enc = preprocessing.OneHotEncoder()
	enc.fit(X)
	NewX =  enc.transform(X).toarray()
	print(NewX[1])
	return NewX,y

def Sklearn(X,y):
	#clf = MLPClassifier(solver='lbfgs', alpha=1e-5,hidden_layer_sizes=(100,),max_iter=100000, random_state=1)
	#clf = svm.SVC(kernel='linear', C=3)
	clf = tree.DecisionTreeClassifier()
	clf.fit(X,y)
	scores = cross_val_score(clf, X, y, cv=5)
	print(scores)
	print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
###############################################
####### main ##############********************
##############################################


#X,y = SingleInput_DAP_Depth("ExampleSequence_Matrix")
#Sklearn(X,y)


X,y = SingleInput("ExampleSequence_Matrix")
Sklearn(X,y)

X,y = MultipleInput("ExampleSequence_Matrix_WindowSize6")
Sklearn(X,y)

X,y = MultipleInput("ExampleSequence_Matrix_WindowSize20")
Sklearn(X,y)
