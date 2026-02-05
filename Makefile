NAME    := glv

CC      := gcc

BASE_CFLAGS := -Wall -Wextra -Werror -std=c11
DEBUG_FLAGS := -g3 -O0 -fsanitize=address,undefined,leak -fno-omit-frame-pointer


CFLAGS  := $(BASE_CFLAGS)
LDFLAGS :=
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
EC_GLV_demo.c \
main.c

OBJS := $(SRCS:.c=.o)

all: $(NAME)

# ---------- RELEASE ----------
$(NAME): $(OBJS)
	$(CC) $(OBJS) -o $(NAME) $(LDFLAGS) $(LIBS)

# ---------- DEBUG ----------
debug: CFLAGS += $(DEBUG_FLAGS)
debug: LDFLAGS += -fsanitize=address,undefined
debug: re

# ---------- RULES ----------
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re debug