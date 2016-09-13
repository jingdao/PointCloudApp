#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <vector>

int num_vocab = 5;
float mse = 0;
bool debug = false;

float getDiff(float* d1, float* d2, int descSize) {
	float f = 0;
	for (int i=0;i<descSize;i++)
		f += (d1[i]-d2[i]) * (d1[i]-d2[i]);
	return f;
}

void getKMeans(float** desc, float** vocab, int* match, int numDesc, int descSize, int k) {
	float** newVocab = new float*[k];
	for (int i=0;i<k;i++) {
		newVocab[i] = new float[descSize];
		memcpy(vocab[i], desc[rand() % numDesc], descSize * sizeof(float));
	}
	int* count = new int[k];
	bool updated = true;
	while (updated) {
		updated = false;
		for (int i=0;i<k;i++)
			memset(newVocab[i],0,descSize*sizeof(float));
		memset(count,0,k*sizeof(int));
		for (int i=0;i<numDesc;i++) {
			int minID=0;
			float minD=getDiff(desc[i],vocab[0],descSize);
			for (int j=1;j<k;j++) {
				float d = getDiff(desc[i],vocab[j],descSize);
				if (d < minD) {
					minD = d;
					minID = j;
				}
			}
			for (int j=0;j<descSize;j++)
				newVocab[minID][j] += desc[i][j];
			count[minID]++;
			if (minID != match[i]) {
				match[i] = minID;
				updated = true;
			}
		}
		for (int i=0;i<k;i++) {
			if (count[i] == 0) {
				memcpy(vocab[i], desc[rand() % numDesc], descSize * sizeof(float));
			} else {
				for (int j=0;j<descSize;j++)
					vocab[i][j] = newVocab[i][j] / count[i];
			}
		}
	}
	for (int i=0;i<k;i++)
		delete[] newVocab[i];
	delete[] newVocab;
	delete[] count;
}

bool loadFileToDictionary(char* name, std::vector<float*> *dict, int* numDesc, int* descSize) {
	FILE* f = fopen(name,"r");
	if (!f)
		return false;
	char buf[2048];
	while (fgets(buf, 2048, f)) {
		if (sscanf(buf, "POINTS %d", numDesc) == 1) {
		} else if (sscanf(buf,"COUNT %d",descSize) == 1) {
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			break;
		}
	}
	for (int i=0;i<*numDesc;i++) {
		fgets(buf,2048,f);
		char* c = buf;
		float* desc = new float[*descSize];
		for (int j=0;j<*descSize;j++)
			desc[j] = strtod(c,&c);
		dict->push_back(desc);
	}
	fclose(f);
	return true;
}

int classify1NN(float* desc, float** bags, int* labels, int numLabels, int descSize) {
	int minID=0;
	float minD = getDiff(desc,bags[0],descSize);
	for (int j=1;j<numLabels;j++) {
		float d = getDiff(desc,bags[j],descSize);
		if (d<minD) {
			minD = d;
			minID = j;
		}
	}
	return labels[minID];
}

int validate1NN(int valID, float* desc, float** bags, int* labels, int numLabels, int descSize) {
	int minID=-1;
	float minD;
	for (int j=0;j<numLabels;j++) {
		if (valID == j)
			continue;
		float d = getDiff(desc,bags[j],descSize);
		if (minID<0 || d<minD) {
			minD = d;
			minID = j;
		}
	}
	return labels[minID];
}

int classifyVoting(float* desc, float** bags, int* labels, int numLabels, int descSize, int* maxCount, int numCategories, float threshold) {
	int majority = 0;
	*maxCount = 0;
	int* count = new int[numCategories]();
	for (int i=0;i<numLabels;i++) {
		float d = getDiff(desc,bags[i],descSize);
		if (d < threshold) {
			count[labels[i]-1]++;
			if (count[labels[i]-1] > *maxCount) {
				*maxCount = count[labels[i]-1];
				majority = labels[i];
			}
		}
	}
	delete[] count;
	return majority;
}

int validateVoting(int valID, float* desc, float** bags, int* labels, int numLabels, int descSize, int* maxCount, int numCategories, float threshold) {
	int majority = 0;
	*maxCount = 0;
	int* count = new int[numCategories]();
	for (int i=0;i<numLabels;i++) {
		if (valID == i)
			continue;
		float d = getDiff(desc,bags[i],descSize);
		if (d < threshold) {
			count[labels[i]-1]++;
			if (count[labels[i]-1] > *maxCount) {
				*maxCount = count[labels[i]-1];
				majority = labels[i];
			}
		}
	}
	delete[] count;
	return majority;
}

int main(int argc, char* argv[] ) {

	if (argc < 4) {
		printf("%s train/ test/ {fpfh,spin} [num_vocab] [MSE] [-v]\n",argv[0]);
		return 1;
	}
	if (argc >= 5)
		num_vocab = atoi(argv[4]);
	if (argc >= 6)
		mse = atof(argv[5]);
	for (int i=4;i<argc;i++)
		if (strncmp(argv[i],"-v",2)==0)
			debug = true;
	srand(0);

	//load train data
	char buf[128];
	std::vector<float*> dictionary;
	std::vector<int> intervals;
	std::vector<int> labels;
	int numCategories=0;
	sprintf(buf,"%s/labels.txt",argv[1]);
	FILE* labelFile = fopen(buf,"r");
	if (!labelFile) {
		printf("Cannot open %s\n",buf);
		return 1;
	}
	while (fgets(buf,128,labelFile)) {
		int l = strtol(buf,NULL,10);
		labels.push_back(l);
		if (l > numCategories)
			numCategories = l;
	}
	fclose(labelFile);
	int id=0;
	int descSize=0;
	while (true) {
		int numDesc=0;
		sprintf(buf,"%s/%d-cloud.pcd-%s.pcd",argv[1],id,argv[3]);
		if (loadFileToDictionary(buf,&dictionary,&numDesc,&descSize)) {
			intervals.push_back(numDesc);
			id++;
		} else
			break;
	}
	if (debug)
		printf("Loaded %lu train descriptors (dimension %d)\n",dictionary.size(),descSize);

	//load test data
	std::vector<float*> test_dictionary;
	std::vector<int> test_intervals;
	std::vector<int> test_labels;
	sprintf(buf,"%s/labels.txt",argv[2]);
	labelFile = fopen(buf,"r");
	if (!labelFile) {
		printf("Cannot open %s\n",buf);
		return 1;
	}
	while (fgets(buf,128,labelFile)) {
		int l = strtol(buf,NULL,10);
		test_labels.push_back(l);
	}
	fclose(labelFile);
	id=0;
	descSize=0;
	while (true) {
		int numDesc=0;
		sprintf(buf,"%s/%d-cloud.pcd-%s.pcd",argv[2],id,argv[3]);
		if (loadFileToDictionary(buf,&test_dictionary,&numDesc,&descSize)) {
			test_intervals.push_back(numDesc);
			id++;
		} else
			break;
	}
	if (debug)
		printf("Loaded %lu test descriptors (dimension %d)\n",test_dictionary.size(),descSize);

	
	//perform K Means to get Vocabulary
	float** vocab = new float*[num_vocab];
	for (int i=0;i<num_vocab;i++)
		vocab[i] = new float[descSize];
	int* match = new int[dictionary.size()]();
	getKMeans(dictionary.data(),vocab,match,dictionary.size(),descSize,num_vocab);

	//compute Bag-of-Words
	float** bags = new float*[labels.size()];
	for (int i=0;i<labels.size();i++)
		bags[i] = new float[num_vocab]();
	for(size_t i=0,n=0;i<intervals.size();i++) {
		int count=0;
		for (int j=0;j<intervals[i];j++,n++) {
			int minID=0;
			float minD=getDiff(dictionary[n],vocab[0],descSize);
			for (int k=1;k<num_vocab;k++) {
				float d = getDiff(dictionary[n],vocab[k],descSize);
				if (d < minD) {
					minD = d;
					minID = k;
				}
			}
			bags[i][minID] += 1;
			count++;
		}
		if (debug)
			printf("Obj %2lu (class %d desc %3d):",i,labels[i],intervals[i]);
		for (int j=0;j<num_vocab;j++) {
			bags[i][j] /= count;
			if (debug)
				printf(" %.3f",bags[i][j]);
		}
		if (debug)
			printf("\n");
	}

	//compute Validation Error
	int numCorrect=0;
	for (size_t i=0;i<labels.size();i++) {
		int neighbors = 1, prediction;
		if (mse > 0)
			prediction = validateVoting(i,bags[i],bags,labels.data(),labels.size(),num_vocab,&neighbors,numCategories,mse * num_vocab);
		else
			prediction = validate1NN(i,bags[i],bags,labels.data(),labels.size(),num_vocab);
		if (labels[i] == prediction)
			numCorrect++;
		if (debug)
			printf("Obj %2lu match %2d neighbors %d %s\n",i,prediction,neighbors,labels[i]==prediction ? "\033[32m/\033[0m" : "\033[31mx\033[0m");
	}
	printf("K=%2d validation accuracy: %2d/%2lu (%.2f)\n",num_vocab,numCorrect,labels.size(),1.0 * numCorrect / labels.size());

	//compute Test Error
	numCorrect=0;
	int* numSamples = new int[numCategories]();
	int* TP = new int[numCategories]();
	float* test_bag = new float[num_vocab];
	for(size_t i=0,n=0;i<test_intervals.size();i++) {
		int count=0;
		memset(test_bag,0,num_vocab*sizeof(float));
		int minID;
		float minD;
		for (int j=0;j<test_intervals[i];j++,n++) {
			minID=0;
			minD=getDiff(test_dictionary[n],vocab[0],descSize);
			for (int k=1;k<num_vocab;k++) {
				float d = getDiff(test_dictionary[n],vocab[k],descSize);
				if (d < minD) {
					minD = d;
					minID = k;
				}
			}
			test_bag[minID] += 1;
			count++;
		}
		for (int j=0;j<num_vocab;j++) {
			test_bag[j] /= count;
		}
		int neighbors = 1,prediction;
		if (mse > 0)
			prediction = classifyVoting(test_bag,bags,labels.data(),labels.size(),num_vocab,&neighbors,numCategories,mse * num_vocab);
		else
			prediction = classify1NN(test_bag,bags,labels.data(),labels.size(),num_vocab);
		numSamples[test_labels[i]-1]++;
		if (prediction == test_labels[i]) {
			TP[test_labels[i]-1]++;
			numCorrect++;
		}
		if (debug)
			printf("Test %2lu match %2d neighbors %2d %s\n",i,prediction,neighbors,test_labels[i]==prediction ? "\033[32m/\033[0m" : "\033[31mx\033[0m");
	}
	printf("K=%2d test accuracy: %2d/%2lu (%.2f)\n",num_vocab,numCorrect,test_labels.size(),1.0 * numCorrect / test_labels.size());
	if (debug)
		for (int i=0;i<numCategories;i++)
			printf("class %d (%2d samples): %2d (%.2f)\n",i+1,numSamples[i],TP[i],1.0 * TP[i] / numSamples[i]);

	delete[] match;
	delete[] test_bag;
	for (int i=0;i<labels.size();i++)
		delete[] bags[i];
	delete[] bags;
	for (int i=0;i<num_vocab;i++)
		delete[] vocab[i];
	delete[] vocab;
	for (size_t i=0;i<dictionary.size();i++)
		delete[] dictionary[i];
	for (size_t i=0;i<test_dictionary.size();i++)
		delete[] test_dictionary[i];
}
