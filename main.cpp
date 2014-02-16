#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <xmmintrin.h>

int scoreGenomes(std::vector<std::pair<std::vector<uint8_t>, int>> genomes, bool sse = false);
std::pair<std::vector<uint8_t>, int> scoreTaxa(std::pair<std::vector<uint8_t>, int> str1, std::pair<std::vector<uint8_t>, int> str2);
std::pair<std::vector<uint8_t>, int> scoreTaxaSSE(std::pair<std::vector<uint8_t>, int> str1, std::pair<std::vector<uint8_t>, int> str2);
void printTaxa(std::pair<std::vector<uint8_t>, int> taxa);

__m128i zero, comp, tempOr, first, second, result;

uint8_t *pZero = reinterpret_cast<uint8_t *>(&zero);
uint8_t *pFirst = reinterpret_cast<uint8_t *>(&first);
uint8_t *pSecond = reinterpret_cast<uint8_t *>(&second);
uint8_t *pResult = reinterpret_cast<uint8_t *>(&result);

int main(){
	int LENGTH = 512;
	int NUM_TAXA = 700;

	std::mt19937 rng(time(NULL));
	std::uniform_int_distribution<> dist(0,4);
	std::vector<std::pair<std::vector<uint8_t>, int>> genomes(NUM_TAXA);
	std::vector<std::pair<std::vector<uint8_t>, int>> layer, clone;

	for(int i = 0; i < 16; ++i){
		pZero[i] = 0;
	}

	std::generate(genomes.begin(), genomes.end(),
		[&](){
			std::vector<uint8_t> tempVec(LENGTH);
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

	auto start = std::chrono::steady_clock::now();
	auto score = scoreGenomes(genomes);
	std::cout<<"Score: "<<score<<" Ran in "<<(std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - start).count())<<" seconds"<<std::endl;
	start = std::chrono::steady_clock::now();
	score = scoreGenomes(genomes, true);
	std::cout<<"SSE Score: "<<score<<" Ran in "<<(std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - start).count())<<" seconds"<<std::endl;
}

int scoreGenomes(std::vector<std::pair<std::vector<uint8_t>, int>> genomes, bool sse){
	std::vector<std::pair<std::vector<uint8_t>, int>> layer;
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

std::pair<std::vector<uint8_t>, int> scoreTaxa(std::pair<std::vector<uint8_t>, int> str1, std::pair<std::vector<uint8_t>, int> str2){
	int score = str1.second + str2.second;
	std::vector<uint8_t> out(str1.first.size());
	for(int i = 0; i < str1.first.size(); ++i){
		out[i] = str1.first[i] & str2.first[i];
		if(!out[i]){
			out[i] = str1.first[i] | str2.first[i];
			++score;
		}
	}
	return std::make_pair(out, score);
}

std::pair<std::vector<uint8_t>, int> scoreTaxaSSE(std::pair<std::vector<uint8_t>, int> str1, std::pair<std::vector<uint8_t>, int> str2){
	int score = str1.second + str2.second;
	std::vector<uint8_t> out;
	for(int i = 0; i < str1.first.size() - 3; i += 16){
		for(int j = 0; j < 16; ++j){
			pFirst[j] = str1.first[i + j];
			pSecond[j] = str2.first[i + j];
		}
		result = _mm_and_si128(first, second);
		for(int j = 0; j < 16; ++j){
			score += !(int)pResult[j];
		}
		tempOr = _mm_or_si128(first, second);
		comp = _mm_cmpeq_epi8(zero, result);
		comp = _mm_and_si128(comp, tempOr);
		result = _mm_or_si128(comp, result);
		out.insert(out.begin(), pResult, pResult + 16);
	}
	return std::make_pair(out, score);
}

void printTaxa(std::pair<std::vector<uint8_t>, int> taxa){
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
