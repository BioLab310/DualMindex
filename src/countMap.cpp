
#include "countMap.h"

void countMapForSAI(vector <struct seed_message> VS_message)
{
	map<int,int> map_vs_message;
	map<int,int> :: iterator map_vs_message_Iter;
	cout << VS_message.size()<< endl;
	for(uint32_t i=0;i<VS_message.size();i++)
	{
		uint32_t x;
		if(VS_message[i].saInterval[1]>=VS_message[i].saInterval[0])
		{
			x=VS_message[i].saInterval[1]-VS_message[i].saInterval[0]+1;
		}
		else
		{
			x=0;
		}
		if(map_vs_message.count(x)>0){
			map_vs_message[x]++;
		}else{
			map_vs_message.insert(pair<int,int>(x,1));
		}
	}

	ofstream output("SAI.txt");

	for(map_vs_message_Iter=map_vs_message.begin();map_vs_message_Iter!=map_vs_message.end();map_vs_message_Iter++)
	{
		output << map_vs_message_Iter->first << " " << map_vs_message_Iter->second<<endl;
	}
}
