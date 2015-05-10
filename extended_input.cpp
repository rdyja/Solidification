
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <jaz/string_add.hpp>
#include "extended_input.hpp"

void clear_comment_and_trim(std::string& line) {
	line.erase(std::find(line.begin(), line.end(), '#'), line.end());	
	line = jaz::remove_space(line);
}

std::pair<std::string, std::string> assign(const std::string& line) {
	
	std::vector<std::string> tokens;
        std::back_insert_iterator<std::vector<std::string> > tokens_ii(tokens);
        
        jaz::split('=', line, tokens_ii);
		
	if(tokens.size() != 2)
		throw std::string("Illegal property");
		
	return std::make_pair(tokens[0], tokens[1]);
}

bool group(std::istream& in, MapConf& conf_group) {	
	std::string line;
	while(in) {
		std::getline(in, line);		
		clear_comment_and_trim(line);
		
		if(line.empty())
			continue;
		
		if(line[0] == '}')
			return true;
		conf_group.insert(assign(line));
	}
	return false;
}

bool read_config_file(const std::string& file_name, 
		MapConf& conf, std::vector<MapConf>& groupConfs) {
	std::ifstream file(file_name.c_str());
	
	std::string line;
	while(std::getline(file, line)) {		
		
		clear_comment_and_trim(line);
		if(line.empty())
			continue;
		
		switch(line[0]) {
			case '{':  {
				MapConf mapGroup;
				if(group(file, mapGroup)) {
					groupConfs.push_back(mapGroup);
					break;
				}				
				return false;				
			}
			case '}':
				return false;
			default:
				conf.insert(assign(line));
				break;
		}			
	}
	
	return true;
}

void printMap::operator()(const MapConf& mp) {
		std::cout << "Group number " << ind << std::endl;
		std::cout << "{\n";
		for(auto p : mp) {
			std::cout << "\t" << p.first << " " << p.second << std::endl;
		}
		std::cout << "}\n";
		++ind;
}
