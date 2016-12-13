SRC_DIR = src/
C_DIGEST = digest.c
SRC_DIGEST = $(addprefix $(SRC_DIR), $(C_DIGEST))

FLAGS = -std=c99 -O3
#FLAGS = -std=c99 -g

all: digest

digest: $(SRC_DIGEST)
	gcc $(FLAGS) $(SRC_DIGEST) -o $@
