#ifndef EXTENDED_INPUT_HPP
#define	EXTENDED_INPUT_HPP

#include <map>

typedef std::map<std::string, std::string> MapConf;

bool read_config_file(const std::string& file_name, 
		MapConf& conf, std::vector<MapConf>& groupConfs);

struct printMap {	
	printMap() : ind(0){}
	void operator()(const MapConf& mp);
private:
	unsigned ind;
};

#endif	/* EXTENDED_INPUT_HPP */

