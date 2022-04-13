APP_NAME=dct
OBJS += dct.o

CXX = g++ -m64 -std=c++1y
CXXFLAGS = -I. -O3 -Wall -fopenmp -Wno-unknown-pragmas

default: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lX11

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@ -lX11

clean:
	rm -rf $(APP_NAME) $(OBJS)
