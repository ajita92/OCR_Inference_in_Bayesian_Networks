#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <cassert>
#include <map>
#include "Factor.cpp"

using namespace std;
double time1 = 0.0;
double LikeLihood = 0;
map <char, int> char2int;
map <int, char> int2char;
vector <string> correctWords;
long totalChars = 0, charsSame = 0; 
long totalWords = 0, wordsSame = 0;
double ocr[1000][10], transition[10][10];
char charaters[10] = {'e', 't', 'a', 'o', 'i', 'n', 's', 'h', 'r', 'd'};

struct node
{
	Factor* psi;
	Factor* tao; // generally tao
	node* parent;

	vector <long> psi_scope, tao_scope; // generally called scope
	vector <node*> children;
	long imgId_toeliminate, id;
};

node *createNode(long _imgIDtoeliminate, vector <long> scope, long scopeID)
{
	node *newNode =  new node;
	newNode->imgId_toeliminate = _imgIDtoeliminate;
	newNode->parent = NULL;
	newNode->psi_scope = scope;
	newNode->id = scopeID;
	newNode->psi = NULL;
	newNode->tao = NULL;
	scope.erase(remove(scope.begin(), scope.end(), scopeID), scope.end());
	newNode->tao_scope = scope;
	return newNode;
}

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
void updateMatrix(long** matrix, vector <long> image, long size, long begin, 
	vector <Factor*> &initial_facts)
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

				matrix[i][j] = 1;
				matrix[j][i] = 1;

				Factor* f = createFactor(vec, 10);
				computeTable(vec.size(), f, tmp, 3, vec, 0);
				initial_facts.push_back(f);

				vec.pop_back();
				vec.pop_back();
			}
		}
	}
}
void updatePairs(long** matrix, vector <long> image, long size1, 
	vector <Factor*> &initial_facts)
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

				matrix[i][j] = 1;
				matrix[j][i] = 1;

				Factor* f = createFactor(vec, 10);
				computeTable(vec.size(), f, tmp, 4, vec, 0);
				initial_facts.push_back(f);

				vec.pop_back();
				vec.pop_back();
			}
		}
	}
}
void updateTans(long** matrix, long size, long begin, vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for (long i = begin; i < size - 1; i++)
	{
		vec.push_back(i);
		vec.push_back(i + 1);

		matrix[i][i + 1] = 1;
		matrix[i + 1][i] = 1;

		Factor* f = createFactor(vec, 10);
		computeTable(vec.size(), f, tmp, 2, vec, 0);
		initial_facts.push_back(f);

		vec.pop_back();
		vec.pop_back();
	}
}

void initialize_matrix(long** matrix, vector <long> image, long size1, long size2, 
	long category, vector <Factor*> &initial_facts)
{
	if (category >= 1)
	{
		vector <long> vec, tmp;
		for (long i = 0; i < image.size(); i++)
		{
			vec.push_back(i);
			Factor *f = createFactor(vec, 10);
			computeTable(vec.size(), f, tmp, 1, vec, image[i]);
			initial_facts.push_back(f);
			vec.pop_back();
		}
	}

	if (category >= 2)
	{
		updateTans(matrix, size1, 0, initial_facts);
		updateTans(matrix, image.size(), size1, initial_facts);
	}

	if (category >= 3)
	{
		updateMatrix(matrix, image, size1, 0, initial_facts);
		updateMatrix(matrix, image, image.size(), size1, initial_facts);
	}

	if (category >= 4)
		updatePairs(matrix, image, size1, initial_facts);

}

bool allmarked(vector <bool> &marked)
{
	for (long m = 0; m < marked.size(); m++)
	{
		if (marked[m] == false) return false;
	}
	return true;
}
long choose_Index(long** matrix, vector <bool> marked)
{
	long minFillIndex = -1;
	long minfillEdges = 1000000;
	vector <long> fillnbrs;

	for (long i = 0; i < marked.size(); i++)
	{
		if (marked[i] == false)
		{
			vector <long> nbrs;
			long filledges = 0;

			for (long j = 0; j < marked.size(); j++)
			{
				if (matrix[i][j] == 1 && marked[j] == false) nbrs.push_back(j);
			}

			for (long n = 0; n < nbrs.size(); n++)
			{
				for (long m = n + 1; m < nbrs.size(); m++)
				{
					if (matrix[nbrs[n]][nbrs[m]] == 0) filledges++;
				}
			}

			if (filledges < minfillEdges) 
			{
				fillnbrs.clear();
				fillnbrs = nbrs;
				minfillEdges = filledges;
				minFillIndex = i;
			}
		}
	}
	for (long f = 0; f < fillnbrs.size(); f++)
	{
		for (long f1 = f + 1; f1 < fillnbrs.size(); f1++)
		{
			matrix[fillnbrs[f]][fillnbrs[f1]] = 1;
			matrix[fillnbrs[f1]][fillnbrs[f]] = 1;
		}
	}
	return minFillIndex;
}

vector <long> greedy_Ordering(long** matrix, long size)
{
	long index;
	vector <bool> marked;
	vector <long> ordering;

	for (long i = 0; i < size; i++)
		marked.push_back(false);

	while(!allmarked(marked))
	{
		index = choose_Index(matrix, marked);
		marked[index] = true;
		ordering.push_back(index);
	}
	return ordering;

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

void printRoot(node* root)
{
	if (root == NULL) return;

	cout << root->id << " : ";
	for (auto a : root->psi_scope)
		cout << a << " ";

	if (root->parent != NULL)
	{
		cout << "||" ;
		cout << root->parent->id << " : ";
		for (auto a : root->parent->psi_scope)
			cout << a << " ";
	}
	cout << endl;
	for (auto a : root->children)
		printRoot(a);
}

void printTree(vector <node*> nodeSet)
{
	for (auto a : nodeSet)
	{
		printRoot(a);
	}
}

bool ispr(long scpID, vector <long> scp)
{
	for (auto a : scp)
	{
		if (scpID == a) return true;
	}
	return false;
}

void init_cliqueTree(node *root, vector <Factor*> &initial_facts)
{
	vector <Factor*> mul_factors;
	for (auto f : initial_facts)
	{
		if (issubset(f->scope, root->psi_scope) && ispr(root->id, f->scope)) mul_factors.push_back(f);
	}

	for (auto i : mul_factors)
	{
		initial_facts.erase(remove(initial_facts.begin(), initial_facts.end(), i), initial_facts.end());
	}

	root->psi = compute_psi(mul_factors, 10);

	//cout << "root psi" << endl;
	//printFactor(root->psi);
	for (auto n : root->children)
	{
		init_cliqueTree(n, initial_facts);
	}
}

void upwardPass(node *root)
{
	for (auto n : root->children)
	{
		upwardPass(n);
		assert(root->psi != NULL);
		root->psi = mutiply(root->psi, n->tao, 10);
	}

	root->tao = marginalize(root->psi, root->id, 10);

}

void downwardPass(node *root)
{
	if (root->parent != NULL)
	{
		Factor *sum_out;
		sum_out = sumOutVariables(root->parent->psi, root->parent->psi->scope, root->tao->scope, 10);
	
		sum_out = Divide(sum_out, root->tao, 10);
		root->psi = mutiply(root->psi, sum_out, 10);
	}
	for (auto r : root->children)
		downwardPass(r);

}

void Belief_Propagation(vector <node*> &nodeSet, vector <Factor*> &initial_facts)
{
	for (auto a : nodeSet)
	{
		init_cliqueTree(a, initial_facts);
		//assert(initial_facts.size() == 0);
		upwardPass(a);
		downwardPass(a);

	}
}

void pred(node* root, vector <long> &vec_str, string finalWord)
{

	vector <long> vec;
	vec.push_back(root->id);
	char ch = finalWord[root->id];
	Factor* fact = sumOutVariables(root->psi, root->psi->scope, vec, 10);

	//printFactor(fact);
	long maxnum = 0, assign = fact->table[0].first[0];
	double num = 0, Z = exp(fact->table[0].second);

	if (fact->table[0].first[0] == char2int[ch])
		num = exp(fact->table[0].second);

	for(long i = 1; i < fact->table.size(); i++)
	{
		if (fact->table[i].first[0] == char2int[ch])
			num = exp(fact->table[i].second);

		Z += exp(fact->table[i].second);
		if (fact->table[i].second > fact->table[maxnum].second)
		{
			maxnum = i;
			assign = fact->table[i].first[0];
		}
	}

	LikeLihood += log((double)num/Z);
	vec_str[root->id] = assign;

	for (auto r : root->children)
		pred(r, vec_str, finalWord);
}

string prediction(vector <node*> &nodeSet, vector <long> &vec_str, string finalWord)
{
	string str = "";
	long start_s = clock();

	for (auto a : nodeSet)
		pred(a, vec_str, finalWord); 

	for (auto i : vec_str)
		str += int2char[i];
	long stop_s = clock();
	time1 += 1.0 * (stop_s - start_s) / CLOCKS_PER_SEC;
	return str;
}

// 1: OCR , 2 : OCR + Trans, 3 : OCR + Trans + Skip, 4 : OCR + Trans + Skip + Parskip
void Modelling(vector <long> image, long size1, long size2, string word1, string word2, long category)
{

	vector <long> ordering;
	vector <node*> clique_tree;
	vector <node*> nodeSet;
	vector <node*> removeSet;
	vector <Factor*> initial_facts;

	long size = image.size();
	long** matrix = new long*[size];

	for (long i  = 0; i < size; i++)
	{
		matrix[i] = new long[size];
		for (long j = 0; j < size; j++)
			matrix[i][j] = 0;
	}

	initialize_matrix(matrix, image, size1, size2, category, initial_facts);
	
	//cout << "initial factors size" << initial_facts.size() << endl;
	ordering = greedy_Ordering(matrix, size); // indexes of imageIDs

	for (long o = 0; o < size; o++)
	{
		vector <long> scp;
		// ORDERING.............
		//cout << image[ordering[o]] << " ";

		for (long i = 0; i < size; i++)
		{
			if (matrix[ordering[o]][i] == 1)
				scp.push_back(i);

			matrix[ordering[o]][i] = 0;
			matrix[i][ordering[o]] = 0;
		}
		scp.push_back(ordering[o]);

		node *newNode = createNode(image[ordering[o]], scp, ordering[o]);

		for (long j = 0; j < nodeSet.size(); j++)
		{
			if (issubset(nodeSet[j]->tao_scope, scp)&& ispr(ordering[o], nodeSet[j]->tao_scope))
			{
				newNode->children.push_back(nodeSet[j]);
				nodeSet[j]->parent = newNode;
				removeSet.push_back(nodeSet[j]);
			}
		}
		nodeSet.push_back(newNode);
		for (auto r : removeSet)
		{
			nodeSet.erase(remove(nodeSet.begin(), nodeSet.end(), r), nodeSet.end());
		}

	}

	//PRINT CLIQUE TREE..
	//printTree(nodeSet);

	Belief_Propagation(nodeSet, initial_facts);

	vector <long> vec_str;
	for (long i = 0; i < image.size(); i++)
		vec_str.push_back(-1);
	string str = prediction(nodeSet, vec_str, word1 + word2);
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
	//PRINTING THE WORDS
	//cout << str1 << " "<< str2 << endl;
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