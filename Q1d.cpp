#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cassert>
#include <map>
#include <utility>
#include "Factor.cpp"

using namespace std;
double time1 = 0;
double LikeLihood = 0;
map <char, int> char2int;
map <int, char> int2char;
vector <string> correctWords;
long totalChars = 0, charsSame = 0; 
long totalWords = 0, wordsSame = 0;
double ocr[1000][10], transition[10][10];
char charaters[10] = {'e', 't', 'a', 'o', 'i', 'n', 's', 'h', 'r', 'd'};

void computeTable(long N, Factor *f, vector <long> list, long task, vector <long> vec, long imgID)
{
	if(N == 0)
	{
		double p = 0;
		if (task == 1) p = ocr[imgID][list[0]];
		else if (task == 2) p = transition[list[0]][list[1]];
		else p = (list[0] == list[1]) ? log(5) : 0; // DANGER

		f->table.push_back(make_pair(list, p));
		return;
	}
	else
	{
		for (long i = 0; i < 10; i++)
		{
			list.push_back(i);
			computeTable(N - 1, f, list, task, vec, imgID);
			list.pop_back();
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void updateMatrix(vector <long> image, long size, long begin, vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for(long i = begin; i < size; i++)
	{
		for (long j = i + 1; j < size; j++)
		{
			if(image[i] == image[j])
			{
				vec.push_back(i);
				vec.push_back(j);

				Factor* f = createFactor(vec, 10);
				computeTable(vec.size(), f, tmp, 3, vec, 0);
				initial_facts.push_back(f);

				vec.pop_back();
				vec.pop_back();
			}
		}
	}
}
void updatePairs(vector <long> image, long size1, vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for(long i = 0; i < size1; i++)
	{
		for (long j = size1; j < image.size(); j++)
		{
			if(image[i] == image[j])
			{
				vec.push_back(i);
				vec.push_back(j);

				Factor* f = createFactor(vec, 10);
				computeTable(vec.size(), f, tmp, 4, vec, 0);
				initial_facts.push_back(f);

				vec.pop_back();
				vec.pop_back();
			}
		}
	}
}
void updateTans(long size, long begin, vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for (long i = begin; i < size - 1; i++)
	{
		vec.push_back(i);
		vec.push_back(i + 1);

		Factor* f = createFactor(vec, 10);
		computeTable(vec.size(), f, tmp, 2, vec, 0);
		initial_facts.push_back(f);

		vec.pop_back();
		vec.pop_back();
	}
}

void initialize_Factors(vector <long> image, long size1, long size2, 
	long category, vector <Factor*> &initial_facts, vector <Factor*> &var_facts)
{
	if (category >= 1)
	{
		vector <long> vec, tmp;
		for (long i = 0; i < image.size(); i++)
		{
			vec.push_back(i);
			Factor *f = createFactor(vec, 10);
			computeTable(vec.size(), f, tmp, 1, vec, image[i]);
			var_facts.push_back(f);
			vec.pop_back();
		}
	}

	if (category >= 2)
	{
		updateTans(size1, 0, initial_facts);
		updateTans(image.size(), size1, initial_facts);
	}

	if (category >= 3)
	{
		updateMatrix(image, size1, 0, initial_facts);
		updateMatrix(image, image.size(), size1, initial_facts);
	}

	if (category >= 4)
		updatePairs(image, size1, initial_facts);

}

bool issubset(vector <long> scp1, vector <long> scp2)
{
	bool flag;
	for (long i = 0; i < scp1.size(); i++)
	{
		flag = false;
		for (long j = 0; j < scp2.size(); j++)
		{
			if (scp1[i] == scp2[j])
			{
				flag = true;
				break;
			}
		}
		if (!flag) return false;
	}
	return true;
}
struct node
{
	Factor *Belief, *Belief_old;
	vector <node*> neighbours;
};

node *createNode(Factor* fact)
{
	node *newNode =  new node;
	newNode->Belief = fact;
	newNode->Belief_old = NULL;
	return newNode;
}

void createClusterGraph(vector <node*> &factorNodes, vector <node*> &varNodes, 
	vector <Factor*> &initial_facts, vector <Factor*> &var_facts)
{
	for (long i = 0; i < initial_facts.size(); i++)
	{
		node* newNode = createNode(initial_facts[i]);
		factorNodes.push_back(newNode);
	}

	for (long i = 0; i < var_facts.size(); i++)
	{
		node* newNode = createNode(var_facts[i]);
		varNodes.push_back(newNode);
	}

	for (long j = 0; j < var_facts.size(); j++)
	{
		for (long i = 0; i < initial_facts.size(); i++)
		{
			if (issubset(var_facts[j]->scope, initial_facts[i]->scope))
			{
				varNodes[j]->neighbours.push_back(factorNodes[i]);
				factorNodes[i]->neighbours.push_back(varNodes[j]);
			}
		}
	}
}

void printClusterGraph(vector <node*> &factorNodes, vector <node*> &varNodes)
{
	cout << "Factor Nodes : " <<endl;
	for (long i = 0; i < factorNodes.size(); i++)
	{
		for (auto a : factorNodes[i]->Belief->scope)
			cout << a << " ";
		cout << " :: ";
		for (long j = 0; j < factorNodes[i]->neighbours.size(); j++)
		{
			Factor *f = factorNodes[i]->neighbours[j]->Belief;
			for (auto i : f->scope)
				cout << i << " ";
			if (j != factorNodes[i]->neighbours.size() - 1)
				cout << "--";
		}
		printFactor(factorNodes[i]->Belief);
		cout << endl;
	}
	cout << "Variable Nodes" << endl;
	for (long i = 0; i < varNodes.size(); i++)
	{
		for (auto a : varNodes[i]->Belief->scope)
			cout << a << " ";
		cout << " :: ";
		for (long j = 0; j < varNodes[i]->neighbours.size(); j++)
		{
			Factor *f = varNodes[i]->neighbours[j]->Belief;
			for (auto i : f->scope)
				cout << i << " ";
			if (j != varNodes[i]->neighbours.size() - 1)
				cout << "--";
		}
		printFactor(varNodes[i]->Belief);
		cout << endl;
	}
}

double abs(double a, double b)
{
	double c = a - b;
	if (c < 0) return -1*c;
	return c;
}


Factor *normalize(Factor *f)
{

	double Z = 0;
	for (long i = 0; i < f->tableSize; i++)
		Z += exp(f->table[i].second);

	if (Z != 0)
	{
		for(long i = 0; i < f->tableSize; i++)
			f->table[i].second = log(exp(f->table[i].second)/Z);
	}

	return f;
}


string computeAccuracy(vector <node*> &varNodes, vector <long> image, string finalword)
{
	string str = "";
	vector <long> vec_str;

	long start_s = clock();
	for (long i = 0; i < image.size(); i++)
		vec_str.push_back(-1);

	for (long i = 0; i < varNodes.size(); i++)
	{
		//printFactor(varNodes[i]->Belief);
		long scope = varNodes[i]->Belief->scope[0];
		double maxVal = varNodes[i]->Belief->table[0].second;

		double prob = 0;
		long chk_Assign = char2int[finalword[scope]];
		long assign = varNodes[i]->Belief->table[0].first[0];
		double Z = exp(varNodes[i]->Belief->table[0].second);

		if (chk_Assign == assign) 
			prob = exp(varNodes[i]->Belief->table[0].second);

		for (long j = 1; j < varNodes[i]->Belief->tableSize; j++)
		{
			if (chk_Assign == varNodes[i]->Belief->table[j].first[0])
				prob = exp(varNodes[i]->Belief->table[j].second);

			if (varNodes[i]->Belief->table[j].second >= maxVal)
			{
				maxVal = varNodes[i]->Belief->table[j].second;
				assign = varNodes[i]->Belief->table[j].first[0];
			}
			Z += exp(varNodes[i]->Belief->table[j].second);
		}

		LikeLihood += log((double)prob/Z);
		//LikeLihood += log(prob);
		vec_str[scope] = assign; 
	}

	for (auto i : vec_str)
		str += int2char[i];
	long stop_s = clock();
	time1 += 1.0 * (stop_s - start_s) / CLOCKS_PER_SEC;
	return str;
}

bool difference(Factor *f1, Factor *f2)
{
	for (long i = 0; i < f1->tableSize; i++)
	{
		if (abs(f1->table[i].second, f2->table[i].second) >= 0.0000001) return true;
	}
	return false;
}
void loopyBPAlgorithm(vector <node*> &factorNodes, vector <node*> &varNodes)
{
	vector <long> vec;
	bool flag = false;
	for (auto a : factorNodes)
		a->Belief_old = copyFactor(a->Belief);

	for (auto a : varNodes)
	{
		a->Belief_old = createFactor(a->Belief->scope, 10);
		initialize_table(a->Belief->scope.size(), a->Belief_old , vec, 10, 0);
	}

	long itr = 1;
	while(!flag && itr <= 50)
	{
		//cout << itr << endl;
		itr++;
		flag = true;

		for (long i = 0; i < varNodes.size(); i++)
		{
			
			vector <long> scope = varNodes[i]->Belief->scope;
			for (long j = 0; j < varNodes[i]->neighbours.size(); j++)
			{
				Factor *f = varNodes[i]->neighbours[j]->Belief_old;
				Factor *fact1 = max_sumOutVariables(f, f->scope, scope, 10);
				Factor *fact2 = varNodes[i]->Belief_old;
				Factor *msg = Divide(fact1, fact2, 10);
				varNodes[i]->Belief = mutiply(varNodes[i]->Belief, msg, 10); 
			}
			
			varNodes[i]->Belief = normalize(varNodes[i]->Belief);
			if (difference(varNodes[i]->Belief, varNodes[i]->Belief_old))
				flag = false;
			varNodes[i]->Belief_old = copyFactor(varNodes[i]->Belief);
			
		}
		
		for (long i = 0; i < factorNodes.size(); i++)
		{
			for (long j = 0; j < factorNodes[i]->neighbours.size(); j++)
			{
				vector <long> scope = factorNodes[i]->neighbours[j]->Belief->scope;
				Factor *fact1 = factorNodes[i]->neighbours[j]->Belief_old;
				Factor *f = factorNodes[i]->Belief_old;
				Factor *fact2 = max_sumOutVariables(f, f->scope, scope, 10);
				Factor *msg = Divide(fact1, fact2, 10);

				factorNodes[i]->Belief = mutiply(factorNodes[i]->Belief, msg, 10);
			}
			factorNodes[i]->Belief = normalize(factorNodes[i]->Belief);

			if (difference(factorNodes[i]->Belief, factorNodes[i]->Belief_old))
				flag = false;
			factorNodes[i]->Belief_old = copyFactor(factorNodes[i]->Belief);
		}
	}
}

void Modelling(vector <long> image, long size1, long size2, string word1, string word2, long category)
{
	vector <Factor*> initial_facts, var_facts;
	vector <node*> factorNodes, varNodes;

	long size = image.size();
	initialize_Factors(image, size1, size2, category, initial_facts, var_facts);
	createClusterGraph(factorNodes, varNodes, initial_facts, var_facts);
	//printClusterGraph(factorNodes, varNodes);
	loopyBPAlgorithm(factorNodes, varNodes);

	
	string str = computeAccuracy(varNodes, image, word1 + word2);
	//cout << str << endl;
	string str1 = str.substr(0, size1);
	string str2 = str.substr(size1, size2);

	if (str1 == word1) wordsSame++;
	if (str2 == word2) wordsSame++;

	for (long i = 0; i < str1.size(); i++)
	{
		if (str1[i] == word1[i]) charsSame++;
	}
	for (long j = 0; j < str2.size(); j++)
	{
		if (str2[j] == word2[j]) charsSame++;
	}

}


void printAccuracy()
{
	double charAcc, wordAcc, LL;
	long total = totalWords;
	charAcc = (double)charsSame/totalChars * 100;
	wordAcc = (double)wordsSame/totalWords * 100;
	LL = (double)LikeLihood/total;

	/*
	cout << endl;
	cout << "Char Wise Accuracy : ";
	cout << (double)charsSame/totalChars * 100 << endl;

	cout << "Word Wise Accurcay : ";
	cout << (double)wordsSame/totalWords * 100 << endl;

	cout << "Average Log LikeLihood : ";
	cout << (double)LikeLihood/total<< endl;
	cout << "Time Taken : ";
	cout << time1 << endl;*/
	cout << endl;
	printf("%lf %lf %lf %lfs\n\n", charAcc, wordAcc, LL, time1);
}


//RUN :  ./a.out category fileindex
int main(int argc, char const *argv[])
{

	char ch, ch1;
	long n, index = 0;
	double prob;
	string str, imgID1, imgID2, space;

	long category = atoi(argv[1]);
	long fileindex = atoi(argv[2]);

	ifstream ocrDat, transDat, dataDat, wordsDat;
	ocrDat.open("./OCRdataset-2/potentials/ocr.dat");
	transDat.open("./OCRdataset-2/potentials/trans.dat");

	
	if (fileindex == 1)
	{
		dataDat.open("./OCRdataset-2/data/data-tree.dat");
		wordsDat.open("./OCRdataset-2/data/truth-tree.dat");
	}
	else if (fileindex == 2)
	{
		dataDat.open("./OCRdataset-2/data/data-loops.dat");
		wordsDat.open("./OCRdataset-2/data/truth-loops.dat");
	}
	else if(fileindex == 3)
	{
		dataDat.open("./OCRdataset-2/data/data-treeWS.dat");
		wordsDat.open("./OCRdataset-2/data/truth-treeWS.dat");
	}
	else if (fileindex == 4)
	{
		dataDat.open("./OCRdataset-2/data/data-loopsWS.dat");
		wordsDat.open("./OCRdataset-2/data/truth-loopsWS.dat");
	}
	//dataDat.open("./OCRdataset-2/data/input.dat");
	//wordsDat.open("./OCRdataset-2/data/word.dat");

	
	for (int i = 0; i < 10; i++)
	{
		int2char[i] = charaters[i];
		char2int[charaters[i]] = i;
	}
	for (long i = 0; i < 10000; i++)
	{
		ocrDat >> n >> ch >> prob;
		ocr[n][char2int[ch]] = log(prob);
	}
	for (long i = 0; i < 100; i++)
	{
		transDat >> ch >> ch1 >> prob;
		transition[char2int[ch]][char2int[ch1]] = log(prob);
	}
	while(getline(wordsDat, str))
	{
		if (str.size() != 0)
		{
			correctWords.push_back(str);
			totalChars += str.size();
		}
	}

	totalWords = correctWords.size();

	for (int i = 0; i < 10; i++)
	{
		int2char[i] = charaters[i];
		char2int[charaters[i]] = i;
	}

	while(getline(dataDat, imgID1))
	{
		vector <long> image1, image2, image;
		stringstream img1(imgID1);

		while(img1 >> n)
		{
			image1.push_back(n);
			image.push_back(n);
		}

		getline(dataDat, imgID1);
		stringstream img2(imgID1);

		while(img2 >> n)
		{
			image2.push_back(n);
			image.push_back(n);
		}
		Modelling(image, image1.size(), image2.size(), correctWords[index], correctWords[index + 1], category);
		getline(dataDat, space);
		index += 2;
	}
	printAccuracy();
}