SRC_DIR = src/
C_DIGEST = digest.c
C_ANALYZE = analyze.c
SRC_DIGEST = $(addprefix $(SRC_DIR), $(C_DIGEST))
SRC_ANALYZE = $(addprefix $(SRC_DIR), $(C_ANALYZE))

FLAGS = -std=c99 -O3
#FLAGS = -std=c99 -g

all: digest analyze

digest: $(SRC_DIGEST)
	gcc $(FLAGS) $(SRC_DIGEST) -o $@

analyze: $(SRC_ANALYZE)
	gcc $(FLAGS) $(SRC_ANALYZE) -o $@
