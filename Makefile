NAME    := ec_crypto

CC      := gcc
CFLAGS  := -Wall -Wextra -Werror -std=c11
INCLUDES:= -I.
LIBS    := -lgmp

SRCS := \
EC_add_affine.c \
EC_add_proj.c \
EC_square_and_multiply_affine.c \
EC_square_and_multiply_proj.c \
EC_struct.c \
main.c

OBJS := $(SRCS:.c=.o)

all: $(NAME)

$(NAME): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(NAME) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re
