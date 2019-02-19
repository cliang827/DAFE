CC = g++
CCOPTION  = I./src std=gnu++11 O3
FLAGS     = $(addprefix -,$(CCOPTION))
OBJ       = Main.o Exec.o Validation.o Udl.o None.o Cprr.o RlRecom.o RlSim.o Contextrr.o ReckNNGraph.o RkGraph.o CorrelationGraph.o Effectiveness.o Type.o Time.o TxtFile.o
OBJ_DIR   = obj
SRC_DIR   = src
BUILD_DIR = bin

udlf: $(OBJ)
	$(CC) $(FLAGS)   $(addprefix $(OBJ_DIR)/,$(OBJ)) -o $(BUILD_DIR)/udlf


#Core
Main.o: $(SRC_DIR)/Core/Main.cpp
	mkdir -p $(OBJ_DIR)
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Core/Main.cpp -o $(OBJ_DIR)/Main.o

Exec.o: $(SRC_DIR)/Core/Exec.cpp

	$(CC) $(FLAGS)  -c $(SRC_DIR)/Core/Exec.cpp -o $(OBJ_DIR)/Exec.o

Validation.o: $(SRC_DIR)/Core/Validation.cpp

	$(CC) $(FLAGS)  -c $(SRC_DIR)/Core/Validation.cpp -o $(OBJ_DIR)/Validation.o

#Methods
Udl.o: $(SRC_DIR)/Methods/Udl.cpp

	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/Udl.cpp -o $(OBJ_DIR)/Udl.o

None.o: $(SRC_DIR)/Methods/None.cpp

	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/None.cpp -o $(OBJ_DIR)/None.o

Cprr.o: $(SRC_DIR)/Methods/Cprr.cpp

	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/Cprr.cpp -o $(OBJ_DIR)/Cprr.o

RlRecom.o: $(SRC_DIR)/Methods/RlRecom.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/RlRecom.cpp -o $(OBJ_DIR)/RlRecom.o

RlSim.o: $(SRC_DIR)/Methods/RlSim.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/RlSim.cpp -o $(OBJ_DIR)/RlSim.o

Contextrr.o: $(SRC_DIR)/Methods/Contextrr.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/Contextrr.cpp -o $(OBJ_DIR)/Contextrr.o

ReckNNGraph.o: $(SRC_DIR)/Methods/ReckNNGraph.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/ReckNNGraph.cpp -o $(OBJ_DIR)/ReckNNGraph.o

RkGraph.o: $(SRC_DIR)/Methods/RkGraph.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/RkGraph.cpp -o $(OBJ_DIR)/RkGraph.o

CorrelationGraph.o: $(SRC_DIR)/Methods/CorrelationGraph.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Methods/CorrelationGraph.cpp -o $(OBJ_DIR)/CorrelationGraph.o

#Evaluation
Effectiveness.o: $(SRC_DIR)/Evaluation/Effectiveness.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Evaluation/Effectiveness.cpp -o $(OBJ_DIR)/Effectiveness.o

#Utils
Type.o: $(SRC_DIR)/Utils/Type.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Utils/Type.cpp -o $(OBJ_DIR)/Type.o

Time.o: $(SRC_DIR)/Utils/Time.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Utils/Time.cpp -o $(OBJ_DIR)/Time.o

TxtFile.o: $(SRC_DIR)/Utils/TxtFile.cpp
	
	$(CC) $(FLAGS)  -c $(SRC_DIR)/Utils/TxtFile.cpp -o $(OBJ_DIR)/TxtFile.o


clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(BUILD_DIR)/udlf
