CFLAGS= -Wall -m64 -o2
CXX= g++
CONCORDE=/opt/concorde
ILOG=/opt/ibm/ILOG/CPLEX_Studio_Community1263
CPPFLAGS= -DIL_STD -I$(CONCORDE) -I/opt/qsopt64 -I$(CONCORDE)/INCLUDE -I$(ILOG)/cplex/include -I$(ILOG)/concert/include -lconcorde
CPLEXLIB= /opt/qsopt64/qsopt.a -L/opt/qsopt64 -L$(ILOG)/cplex/lib/x86-64_linux/static_pic -lilocplex -lqsopt -lcplex -L$(ILOG)/concert/lib/x86-64_linux/static_pic -lconcert -L$(CONCORDE) -lconcorde -lm -lpthread

bruno:
		$(CXX) $(CFLAGS) $(CPPFLAGS) src/main.cpp -o bin/main $(CPLEXLIB)

stack:
		$(CXX) src/main.cpp -o bin/main -I/opt/qsopt64 -I$(CONCORDE) -L$(ILOG)/cplex/lib/x86-64_linux/static_pic -L$(ILOG)/concert/lib/x86-64_linux/static_pic -lilocplex -lconcert -lcplex /opt/concorde/concorde.a /opt/qsopt64/qsopt.a -lm -lpthread

all:
		g++ -std=c++11 -I/home/rodrigo/Developments/concorde -w src/main.cpp /home/rodrigo/Developments/concorde/concorde.a ./include/qsopt.a -o bin/main -pthread
		@echo "\n--\nCompilação finalizada...\ndigite 'make run' para executar\n--"
run:
		./bin/main