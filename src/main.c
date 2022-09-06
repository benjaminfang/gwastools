#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "config.h"
#include "version.h"

#include "help.h"
#include "sort.h"
#include "merge.h"


int
main (int argc, char * argv[])
{
    if (argc == 1) {
        printf("besdtool, tools to handle besd files\n"
               "benjaminfang.ol@outlook.com\n"
               "Usage: besdtool [--help][--version] subcommand [--help][options] <arguments>\n"
               "besdtool --help for help infomation\n");
        exit(0);
    }

    if (strcmp(argv[1], "--version") == 0) {
        printf("besdtool version: %s\n", BESDTOOL_VERSION);
        exit(0);
    }

    if (strcmp(argv[1], "--help") == 0) {
        printf("besdtool [--help][--version] subcommand [--help][options] <arguments>\n"
               "--version: print besdtool version information and exit.\n"
               "--help: print this help information and exit.\n"
               "subcommand:\n"
               "help           : print help information.\n"
               "sort           : sort besd file.\n"
               "merge          : merge besd file.\n");
        exit(0);
    }
 
    int subcommand_count = 0;

    subcommand_count += help(argc, argv);
    subcommand_count += sort(argc, argv);
    subcommand_count += merge(argc, argv);

    if (subcommand_count == 0) {
        printf("your subcommand was not recgnized.\n");
        help(argc, argv);
    } else if (subcommand_count > 1) {
        printf("This can no be possible.\n");
    }

    return 0;
}