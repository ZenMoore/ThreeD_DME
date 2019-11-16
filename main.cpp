#include "ThreeD_DME.h"
using namespace std;

void parse_argument(int argc, char *argv[], 
string &inputSinksFileName, ZstTree::TopologyMode &topoMode, string &outputFileName) 
{
    for (int i = 1; i < argc; i++) 
    {
        if (strcmp("-i", argv[i]) == 0) {
            i++;
            inputSinksFileName = argv[i];
        } else if (strcmp("-o", argv[i]) == 0) {
            i++;
            outputFileName = argv[i] ;
        } else if (strcmp("-M", argv[i]) == 0)  { // for topology model
            i++;
            int model;
            sscanf(argv[i], "%d", &model);
            topoMode = (ZstTree::TopologyMode) model;
        } else {
            printf("Argument %d incorrect\n", i);
            exit(0);
        }
    }
}

const string ZSTUSAGE = "usage: zst -i inputFileName -M topoMode -o outputFileName \n"; 

int main(int argc, char *argv[]) 
{
    string inputSinksFileName = "" ; // mandatory
    string outputFileName = ""; // mandatory
    ZstTree::TopologyMode topoMode = ZstTree::GREEDYMODE;

    parse_argument (argc, argv, inputSinksFileName, topoMode, outputFileName);

    string fileNames[]{
        // "../Benchmark/point4_p2.txt",
        // "../Benchmark/point8_p2.txt",
        // "../Benchmark/point16_p2.txt",
        // "../Benchmark/point32_p2.txt",
        // "../Benchmark/point64_p2.txt",
        // "../Benchmark/point128_p2.txt",
        "../Benchmark/r1_p2.txt",
        "../Benchmark/r2_p2.txt",
        "../Benchmark/r3_p2.txt",
        "../Benchmark/r4_p2.txt",
        "../Benchmark/r5_p2.txt"
    };

    int n = sizeof(fileNames) / sizeof(string);

    FILE *f = fopen(outputFileName.c_str(), "w");
    for (int i = 0; i < n; i++) {
        inputSinksFileName = fileNames[i];
        if (inputSinksFileName.empty() || outputFileName.empty()) { // mandatory 
            cout << ZSTUSAGE << endl;
            return -1;
        } else {
            fprintf(f, "%s:\n", inputSinksFileName.c_str());
            DOUBLE totalLen[2][2];
            for (ZstTree::TopologyMode enumi = ZstTree::GREEDYMODE; enumi <= ZstTree::MMMMODE; enumi = (ZstTree::TopologyMode)(enumi + 1)) {
                if (enumi == ZstTree::GREEDYMODE) {
                    fprintf(f, "\tMode: GREEDYMODE\n");
                } else {
                    fprintf(f, "\tMode: MMMMODE\n");
                }
                for (YESNO enumj = NO; enumj <= YES; enumj = (YESNO)(enumj + 1)) {
                    if (enumi == ZstTree::MMMMODE && enumj == YES) {
                        continue;
                    }
                    fprintf(f, "\t\tLOOK_AHEAD: %d\n", enumj);
                    ThreeD_DME dme (inputSinksFileName, enumi, enumj);
                    dme.ConstructTree();
                    DOUBLE temp = dme.TotalLength();
                    fprintf(f, "\t\t\tTotalLength: %lf\n", temp.value);
                    fprintf(f, "\t\t\tTotalNumTsv: %d\n", dme.TotalNumTsv());
                    fprintf(f, "\t\t\tTotalPower:  %lf watt\n", dme.TotalPower().value);
                    totalLen[enumi][enumj] = temp;
                }
            }
            DOUBLE reduction1 = ((totalLen[ZstTree::MMMMODE][0] - totalLen[ZstTree::GREEDYMODE][0]) / totalLen[ZstTree::MMMMODE][0]) * 100;
            DOUBLE reduction2 = ((totalLen[ZstTree::GREEDYMODE][0] - totalLen[ZstTree::GREEDYMODE][1]) / totalLen[ZstTree::GREEDYMODE][0]) * 100;
            fprintf(f, "Reduction1: %.2lf%%\n", reduction1.value);
            fprintf(f, "Reduction2: %.2lf%%\n", reduction2.value);
            fprintf(f, "\n");
        }
    }
    fclose(f);
    return 0;
}
