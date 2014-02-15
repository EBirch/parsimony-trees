#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>

std::pair<std::vector<int>, int> scoreTaxa(std::pair<std::vector<int>, int> str1, std::pair<std::vector<int>, int> str2);
void printTaxa(std::vector<int> taxa);

int main(){
	int LENGTH = 1000;
	int NUM_TAXA = 20;

	std::mt19937 rng(time(NULL));
	std::uniform_int_distribution<> dist(0,4);
	std::vector<std::pair<std::vector<int>, int>> genomes(NUM_TAXA);
	std::vector<std::pair<std::vector<int>, int>> layer;

	std::generate(genomes.begin(), genomes.end(),
		[&](){
			std::vector<int> tempVec(LENGTH);
			std::generate(tempVec.begin(), tempVec.end(),
				[&](){
					int tempHex;
					switch(dist(rng)){
						case 0: tempHex = 0x1;
							break;
						case 1: tempHex = 0x2;
							break;
						case 2:	tempHex = 0x4;
							break;
						case 3: tempHex = 0x8;
							break;
						case 4: tempHex = 0x10;
							break;
					}
					return tempHex;
			});
			return std::make_pair(tempVec, 0);
		});	

	while(genomes.size() > 1){
		layer = genomes;
		genomes.clear();
		// for(int i = 0; i < layer.size(); i += 2){
		while(layer.size() > 0){
			if(layer.size() == 1){
				genomes.push_back(layer[0]);
				continue;
			}
			uint minScore = -1;
			int index = 0;
			for(int i = 1; i < layer.size(); ++i){
				auto temp = scoreTaxa(layer[0], layer[i]);
				if(temp.second < minScore){
					minScore = temp.second;
					index = i;
				}
			}
			genomes.push_back(scoreTaxa(layer[0], layer[index]));
			layer.erase(layer.begin() + index);
			layer.erase(layer.begin());
		}
	}

	std::cout<<"Score: "<<genomes[0].second<<std::endl;
	// for(int i = 0; i < genomes.size() - 1; i += 2){
	// 	auto temp = scoreTaxa(genomes[i], genomes[i + 1]);
	// 	// std::cout<<temp<<std::endl;
	// }
}

std::pair<std::vector<int>, int> scoreTaxa(std::pair<std::vector<int>, int> str1, std::pair<std::vector<int>, int> str2){
	int score = str1.second + str2.second;
	// printTaxa(str1);
	// printTaxa(str2);
	std::vector<int> out(str1.first.size());
	for(int i = 0; i < str1.first.size(); ++i){
		out[i] = str1.first[i] & str2.first[i];
		if(!out[i]){
			out[i] = str1.first[i] | str2.first[i];
			++score;
		}
	}
	// printTaxa(out);
	return std::make_pair(out, score);
}

void printTaxa(std::vector<int> taxa){
	for(auto &gene : taxa){
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
