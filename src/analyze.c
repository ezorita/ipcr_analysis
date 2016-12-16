#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <search.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <ctype.h>

#define HASH_SIZE 2048
#define MAX_FRAGMENT_SIZE 2000

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#define lowercase(s) for(char * p = s;*p;++p) *p=tolower(*p)

// Structures.
typedef struct {
   char * chr;
   int    cnt;
   int  * re_site;
} isdchr_t;

typedef struct {
   char * re_seq;
   int    cut_fw;
   int    cut_rv;
} re_t;

// Function headers.
re_t     parse_isd    (int re_fd, char* re_name, struct hsearch_data* htable);
int      bisection    (int* data, int beg, int end, int target);

int main(int argc, char *argv[])
{
   if (argc != 5) {
      fprintf(stderr, "usage: %s <integration site table> <organism name> <RE1,RE2,RE3...> <construct length>\n", argv[0]);
      exit(1);
   }

   char * organism;
   char * in_file;
   char * re_name_;
   int    cons_len;

   in_file  = argv[1];
   organism = argv[2];
   re_name_ = argv[3];
   cons_len = atoi(argv[4]);


   // Open input file.
   FILE * infile = fopen(in_file, "r");
   if (infile == NULL) {
      fprintf(stderr,"error opening input file: %s.\n", in_file);
      exit(1);
   }

   // Open digest files.
   int re_cnt = 0;
   int re_fd[5];
   char * re_name[5];
   char * re_nm = strtok(re_name_, ",");
   do {
      // To lowercase.
      char * name = strdup(re_nm);
      lowercase(name);
      char * db_path = malloc(strlen(name)+strlen(organism)+9);
      sprintf(db_path, "db/%s/%s.isd", organism, name);

      // Open file and store fd.
      re_name[re_cnt] = name;
      re_fd[re_cnt] = open(db_path, O_RDONLY);
      free(db_path);
      if (re_fd[re_cnt] < 0) {
         fprintf(stderr, "error while opening: %s.\n",db_path);
         exit(1);
      }
      re_cnt++;
   } while ((re_nm = strtok(NULL, ",")));


   // Create a different hash table for each Digestion? Or only one hash table and assign the keys as "re_name+chr_name?". May be easier given the implementation of the hash table in C.
   struct hsearch_data htable = {0};
   hcreate_r(HASH_SIZE, &htable);

   re_t * re_info = malloc(re_cnt*sizeof(re_t));

   // Parse digest files.
   for (int i = 0; i < re_cnt; i++) {
      re_info[i] = parse_isd(re_fd[i], re_name[i], &htable);
   }

   // Read input file (integrations) and perform 're_cnt' bisection searches.
   size_t n = 100;
   char * line = malloc(n);
   ssize_t bytes;
   
   int lineno = 0;
   while ((bytes = getline(&line, &n, infile)) > 0) {
      lineno++;
      // Parse input line.
      char * sqname = strtok(line,"\t\n");
      char *    chr = strtok(NULL,"\t\n");
      int     locus = atoi(strtok(NULL,"\t\n"));
      char * strand = strtok(NULL,"\t\n");

      // Vars.
      int next_line = 0;
      int site = -1;
      int re = -1;
      
      for (int i = 0; i < re_cnt; i++) {
         // Generate hash key.
         char * hkey = malloc(strlen(chr)+strlen(re_name[i])+1);
         sprintf(hkey, "%s%s", re_name[i], chr);
         // Find chromosome in hash table.
         ENTRY * item;
         hsearch_r((ENTRY){.key = hkey}, FIND, &item, &htable);
         if (item == NULL) {
            fprintf(stderr, "warning (line %d): chromosome not found in digestion file: %s. (Ignoring entry)\n", lineno, chr);
            next_line = 1;
            break;
         }
         isdchr_t * ref = (isdchr_t *) item->data;
         // Bisection algorithm.
         int idx = bisection(ref->re_site, 0, ref->cnt-1, locus);
         if (strand[0] == '+') {
            while (locus-ref->re_site[idx] <= strlen(re_info[i].re_seq)) idx--;
            if (re == -1 || ref->re_site[idx] > site) {
               site = ref->re_site[idx];
               re = i;
            }
         } else {
            while (ref->re_site[idx+1]-locus <= strlen(re_info[i].re_seq)) idx++;
            if (re == -1 || ref->re_site[idx+1] < site) {
               site = ref->re_site[idx+1];
               re = i;
            }
         }
      }
      if (next_line) continue;
      // Output result.
      int mol_len = cons_len + (strand[0] == '+' ? locus - site : site - locus);
      fprintf(stdout, "%s\t%s\t%d\t%s\t%s\t%d\t%d\n", sqname, chr, locus, strand, re_name[re], site, mol_len);

   }
      // Output result.



   return 0;
}


int
bisection
(
 int * data,
 int   beg,
 int   end,
 int   target
 )
{
   if (end - beg < 2) return beg;
   int mid = (beg+end)/2;
   if (target < data[mid]) end = mid;
   else if (target > data[mid]) beg = mid;
   else return mid;
   
   return bisection(data,beg,end,target);
}


re_t
parse_isd
(
 int         re_fd,
 char      * re_name,
 struct hsearch_data * htable
)
{
   struct stat sb;
   char * p, * pmap;
   
   // mmap file.
   if (fstat(re_fd, &sb) == -1) {
      fprintf(stderr, "error reading digestion file (fstat).\n");
      exit(1);
   }

   pmap = mmap(0, sb.st_size, PROT_READ, MAP_SHARED, re_fd, 0);
   if (pmap == NULL) {
      fprintf(stderr, "error reading digestion file (mmap).\n");
      exit(1);
   }

   p = pmap;
   // Read RE information.
   re_t reinfo;
   reinfo.re_seq = p;
   p += strlen(reinfo.re_seq)+1;
   reinfo.cut_fw = *((int *)p);
   p += sizeof(int);
   reinfo.cut_rv = *((int *)p);
   p += sizeof(int);

   // Get number of chromosomes.
   int nchrom = *((int *)p);
   p += sizeof(int);
   
   // Parse each chromosome.
   int re_name_len = strlen(re_name);
   for (int i = 0; i < nchrom; i++) {
      // Create isd chromosome entry.
      isdchr_t * isdchr = malloc(sizeof(isdchr_t));

      // Generate hash table key. (REnamechromosome)
      int chrname_len = strlen(p);
      char * hkey = malloc(re_name_len+chrname_len+1);
      sprintf(hkey, "%s%s", re_name, p);
      
      // Insert isdchr in hash table.
      ENTRY * item;
      hsearch_r((ENTRY){.key = hkey, .data = isdchr}, ENTER, &item, htable);
     
      // Chromosome name.
      isdchr->chr = p;
      p += chrname_len+1;
      // Number of RE sites.
      isdchr->cnt = *((int *)p);
      p += sizeof(int);
      // RE site list.
      isdchr->re_site = ((int *)p);
      p += isdchr->cnt*sizeof(int);
   }

   return reinfo;
}
      
