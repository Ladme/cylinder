cylinder: main.c
	gcc main.c -I$(groan) -L$(groan) -D_POSIX_C_SOURCE=200809L -o cylinder -lgroan -lm -std=c99 -pedantic -Wall -Wextra -O3 -march=native

install: cylinder
	cp cylinder ${HOME}/.local/bin
