INT_DIR := integrator

SRC_DIR := $(INT_DIR)/src
INC_DIR := $(INT_DIR)/include
OBJ_DIR := $(INT_DIR)/build
BIN_DIR := .

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)  

CODE = tidymess.exe

CC     := g++
CFLAGS := -Ofast

.PHONY: all clean

all: $(CODE)

tidymess.exe: $(OBJ) | $(BIN_DIR)
	$(CC) $(CFLAGS) $^ -o $@ $(INC_DIR:%=-I%)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@ $(INC_DIR:%=-I%)
    
$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@
    
clean:
	@$(RM) -rv $(OBJ_DIR) $(CODE)
    
-include $(OBJ:.o=.d)

