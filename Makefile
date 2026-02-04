NAME    := glv

CC      := gcc

CFLAGS  := -Wall -Wextra -Werror -std=c11
DEBUG_FLAGS := -g3 -O0 -fsanitize=address,undefined -fno-omit-frame-pointer

INCLUDES:= -I.
LIBS    := -lgmp

SRCS := \
EC_add_affine.c \
EC_add_proj.c \
EC_square_and_multiply_affine.c \
EC_square_and_multiply_proj.c \
double_scalar_multiplication.c \
precompute_table.c \
EC_struct.c \
EC_endo_phi_GLV.c \
EC_GLV.c \
short_vectors.c \
glv_decompose.c \
quadratic_solver.c \
glv_curves.c \
glv_acceleration.c \
EC_DH.c \
EC_GLV_mpz_t.c \
main.c

OBJS := $(SRCS:.c=.o)

all: $(NAME)

# release build
$(NAME): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(NAME) $(LIBS)

# debug build avec sanitizer
debug: CFLAGS += $(DEBUG_FLAGS)
debug: re

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re debug
