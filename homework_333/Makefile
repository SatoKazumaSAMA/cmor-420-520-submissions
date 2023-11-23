CC = g++
CPPFLAGS = -I./include
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin

INC_FILES := $(wildcard $(INC_DIR)/*.hpp)
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

.PHONY: all clean main

all: main

main: $(OBJ_FILES) main.cpp
	mkdir -p $(BIN_DIR)
	$(CC) $^ $(CPPFLAGS) -o $(BIN_DIR)/main

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_FILES)
	mkdir -p $(OBJ_DIR)
	$(CC) -c $< -o $@ $(CPPFLAGS)

clean:
	rm -f $(BIN_DIR)/main $(OBJ_DIR)/*.o *~

