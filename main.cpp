#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <xmmintrin.h>

int scoreGenomes(std::vector<std::pair<std::vector<unsigned int>, int>> genomes, bool sse = false);
std::pair<std::vector<unsigned int>, int> scoreTaxa(std::pair<std::vector<unsigned int>, int> str1, std::pair<std::vector<unsigned int>, int> str2);
std::pair<std::vector<unsigned int>, int> scoreTaxaSSE(std::pair<std::vector<unsigned int>, int> str1, std::pair<std::vector<unsigned int>, int> str2);
void printTaxa(std::pair<std::vector<unsigned int>, int> taxa);

__m128 zero, comp, tempOr, first, second, result;

unsigned int *pZero = reinterpret_cast<unsigned int *>(&zero);
unsigned int *pFirst = reinterpret_cast<unsigned int *>(&first);
unsigned int *pSecond = reinterpret_cast<unsigned int *>(&second);
unsigned int *pResult = reinterpret_cast<unsigned int *>(&result);

int main(){
	int LENGTH = 16;
	int NUM_TAXA = 2;

	std::mt19937 rng(time(NULL));
	std::uniform_int_distribution<> dist(0,4);
	std::vector<std::pair<std::vector<unsigned int>, int>> genomes(NUM_TAXA);
	std::vector<std::pair<std::vector<unsigned int>, int>> layer, clone;

	// std::cout<<sizeof(unsigned int)<<std::endl;
	// std::cout<<sizeof(float)<<std::endl;
	// return 0;

	// int *pZero = reinterpret_cast<int *>(&zero);
	for(int i = 0; i < 4; ++i){
		pZero[i] = 0;
	}

	std::generate(genomes.begin(), genomes.end(),
		[&](){
			std::vector<unsigned int> tempVec(LENGTH);
			std::generate(tempVec.begin(), tempVec.end(),
				[&](){
					int tempHex;
					switch(dist(rng)){
						case 0: tempHex = 0x1; break;
						case 1: tempHex = 0x2; break;
						case 2:	tempHex = 0x4; break;
						case 3: tempHex = 0x8; break;
						case 4: tempHex = 0x10; break;
					}
					return tempHex;
			});
			return std::make_pair(tempVec, 0);
		});

	std::cout<<"Genome 1: ";
	printTaxa(genomes[0]);
	std::cout<<"Genome 2: ";
	printTaxa(genomes[1]);
	auto score = scoreTaxa(genomes[0], genomes[1]);
	auto scoreSSE = scoreTaxaSSE(genomes[0], genomes[1]);
	std::cout<<"Score: ";
	printTaxa(score);
	std::cout<<"Length: "<<score.first.size()<<std::endl;
	std::cout<<score.second<<std::endl;
	std::cout<<"SSE Score: ";
	printTaxa(scoreSSE);
	std::cout<<"Length: "<<scoreSSE.first.size()<<std::endl;
	std::cout<<scoreSSE.second<<std::endl;

	std::cout<<std::endl<<std::endl;
	for(auto &val : score.first){
		std::cout<<val<<" ";
	}
	std::cout<<std::endl<<std::endl;
	for(auto &val : scoreSSE.first){
		std::cout<<val<<" ";
	}
	std::cout<<std::endl<<std::endl;
	// auto start = std::chrono::steady_clock::now();
	// auto score = scoreGenomes(genomes);
	// std::cout<<"Score: "<<score<<" Ran in "<<(std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - start).count())<<" seconds"<<std::endl;
	// start = std::chrono::steady_clock::now();
	// score = scoreGenomes(genomes, true);
	// std::cout<<"SSE Score: "<<score<<" Ran in "<<(std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - start).count())<<" seconds"<<std::endl;
}

int scoreGenomes(std::vector<std::pair<std::vector<unsigned int>, int>> genomes, bool sse){
	std::vector<std::pair<std::vector<unsigned int>, int>> layer;
	while(genomes.size() > 1){
		layer = genomes;
		genomes.clear();
		while(layer.size() > 0){
			if(layer.size() == 1){
				genomes.push_back(layer[0]);
				continue;
			}
			uint minScore = -1;
			int index = 0;
			for(int i = 1; i < layer.size(); ++i){
				auto temp = sse ? scoreTaxaSSE(layer[0], layer[i]) : scoreTaxa(layer[0], layer[i]);
				if(temp.second < minScore){
					minScore = temp.second;
					index = i;
				}
			}
			genomes.push_back(sse ? scoreTaxaSSE(layer[0], layer[index]) : scoreTaxa(layer[0], layer[index]));
			layer.erase(layer.begin() + index);
			layer.erase(layer.begin());
		}
	}
	return genomes[0].second;
}

std::pair<std::vector<unsigned int>, int> scoreTaxa(std::pair<std::vector<unsigned int>, int> str1, std::pair<std::vector<unsigned int>, int> str2){
	int score = str1.second + str2.second;
	std::vector<unsigned int> out(str1.first.size());
	for(int i = 0; i < str1.first.size(); ++i){
		out[i] = str1.first[i] & str2.first[i];
		if(!out[i]){
			out[i] = str1.first[i] | str2.first[i];
			++score;
		}
	}
	return std::make_pair(out, score);
}

std::pair<std::vector<unsigned int>, int> scoreTaxaSSE(std::pair<std::vector<unsigned int>, int> str1, std::pair<std::vector<unsigned int>, int> str2){
	int score = str1.second + str2.second;
	// __m128 first, second, result;
	// unsigned int *pFirst = reinterpret_cast<unsigned int *>(&first);
	// unsigned int *pSecond = reinterpret_cast<unsigned int *>(&second);
	// unsigned int *pResult = reinterpret_cast<unsigned int *>(&result);
	// std::vector<unsigned int> out(str1.first.size());
	std::vector<unsigned int> out;

	for(int i = 0; i < str1.first.size() - 3; i += 4){
		for(int j = 0; j < 4; ++j){
			pFirst[j] = str1.first[i + j];
			pSecond[j] = str2.first[i + j];
		}
		result = _mm_and_ps(first, second);
		for(int j = 0; j < 4; ++j){
			score += !(int)pResult[j];
		}
		tempOr = _mm_or_ps(first, second);
		//comp = _mm_compeq(zero, result);
		comp = _mm_cmpeq_ps(zero, result);
		comp = _mm_and_ps(comp, tempOr);
		result = _mm_or_ps(comp, result);
		out.insert(out.begin(), pResult, pResult + 4);
	}

	return std::make_pair(out, score);
}

void printTaxa(std::pair<std::vector<unsigned int>, int> taxa){
	for(auto &gene : taxa.first){
		switch(gene){
			case 0x1 : std::cout<<'A'; break;
			case 0x2 : std::cout<<'C'; break;
			case 0x4 : std::cout<<'G'; break;
			case 0x8 : std::cout<<'T'; break;
			case 0x10 : std::cout<<'-'; break;
			default : std::cout<<'[';
				for(int i = 1; i <= 16; i *= 2){
					switch(gene & i){
						case 0x1 : std::cout<<'A'; break;
						case 0x2 : std::cout<<'C'; break;
						case 0x4 : std::cout<<'G'; break;
						case 0x8 : std::cout<<'T'; break;
						case 0x10 : std::cout<<'-'; break;
					}
				}
				std::cout<<']';
		}
	}
	std::cout<<std::endl;
}
