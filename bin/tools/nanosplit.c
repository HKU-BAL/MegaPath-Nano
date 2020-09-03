#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <cstring>
#include <sstream>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

void show_help_and_die(char* my_name)
{
	fprintf(stderr, "Usage: %s <read list> <fastq/fasta file1> <fastq/fasta file2> ...\n", my_name);
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc < 3) show_help_and_die(argv[0]);
	map<string, set<string>> read2file;
	map<string, int> file2index;
	//map<string, set<string>> file2readLists;
	FILE * readList = fopen(argv[1], "r");
	string line;
	size_t len;
	char *ptr = NULL;
	istringstream iss;
	string readId, fileName;
	int numOfDestFile = 0;
	while(getline(&ptr, &len, readList) != -1)
	{
		line = ptr;
		iss.str(line);
		iss >> readId >> fileName;
		read2file[readId].insert(fileName);
		if(file2index.find(fileName) == file2index.end())
		{
			file2index[fileName]=numOfDestFile++;
		}
		//file2readLists[fileName].insert(readId);
	}
	fclose(readList);

	FILE **out;
	out = new FILE*[numOfDestFile];
	char outputFileName[256];
	for(map<string, int>::iterator it = file2index.begin(); it != file2index.end(); it++)
	{
		strcpy(outputFileName, it->first.c_str());
		out[it->second]=fopen(outputFileName, "w");
	}

	gzFile fp;
	kseq_t *seq;
	for(int i=2;i<argc;i++)
	{
		//FILE *out = fopen(outputFileName, "w");
		fp = gzopen(argv[i], "r");
		seq = kseq_init(fp);
		while(kseq_read(seq) >= 0)
		{
			if(read2file.find(seq->name.s) != read2file.end())
			{
				for(set<string>::iterator it = read2file[seq->name.s].begin(); it != read2file[seq->name.s].end(); it++)
				{
					fprintf(out[file2index[*it]], "%c", seq->is_fastq ? '@' : '>');
					fprintf(out[file2index[*it]], "%s", seq->name.s);
					if(seq->comment.l)
						fprintf(out[file2index[*it]], " %s", seq->comment.s);
					fprintf(out[file2index[*it]], "\n");
					fprintf(out[file2index[*it]], "%s\n", seq->seq.s);
					if(seq->is_fastq)
					{
						fprintf(out[file2index[*it]], "+\n");
						fprintf(out[file2index[*it]], "%s\n", seq->qual.s);
					}
				}
			}
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	for(map<string, int>::iterator it = file2index.begin(); it != file2index.end(); it++)
	{
		fclose(out[it->second]);
	}
	return 0;
}