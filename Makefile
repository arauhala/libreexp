
#GPP          = /usr/lib/gcc-snapshot/bin/g++
GPP          = ccache clang++ -Qunused-arguments

CXX 		 =  $(GPP)
CXXLINKER    = $(GPP)

LINKERFLAGS  = 
CXXFLAGS 	 = -Wall -Werror -MMD -MP -fmessage-length=0 -Isrc -std=c++11

SRCS = $(shell find src -name *.cpp)

HEADERS = $(shell find src -name *.h)

TEST_SRCS = $(shell find test -name *.cpp)

GCHS =      

OBJS =		$(patsubst %.cpp,%.o,$(SRCS))

DEPS = 	    $(patsubst %.o,%.d,$(OBJS))

GCHS_DEPS = $(patsubst %.h.gch,%.h.d,$(GCHS))

TEST_OBJS = $(patsubst %.cpp,%.o,$(TEST_SRCS))

TEST_DEPS = $(patsubst %.o,%.d,$(TEST_OBJS))

LIBS =

TARGET =	  libreexp.a

OPTFLAGS = -Ofast -march=native -DNDEBUG

TEST_TARGET = libreexptest

def : debug

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += $(OPTFLAGS)
release: all

fastbuild: all

optdebug: CXXFLAGS += $(OPTFLAGS) -g
optdebug: all

clean_profile: 
	rm -f gmon.out

profile: CXXFLAGS += -pg
profile: LINKERFLAGS += -pg
profile: all

testrun: 
	./libexptest

libexptest.profile : clean_profile profile testrun
	gprof ./libexptest > libexptest.profile

leafpad_profile: libexptest.profile
	leafpad libexptest.profile &

optprofile: CXXFLAGS += $(OPTFLAGS) -pg
optprofile: LINKERFLAGS += -pg
optprofile: all

libexptest.optprofile: clean_profile optprofile testrun
	gprof ./libexptest > libexptest.optprofile

leafpad_optprofile: libexptest.optprofile
	leafpad libexptest.optprofile &

$(GCHS): %.h.gch: %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@  

$(OBJS): %.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

$(TEST_OBJS): %.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -Itest $< -o $@
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

-include $(DEPS) $(TEST_DEPS) $(GCHS_DEPS) 

$(TARGET):	$(OBJS) $(GCHS)
	$(AR) rcs $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TEST_OBJS) $(DEPS) $(TEST_DEPS) $(TARGET) $(TEST_TARGET) $(GCHS) $(GCHS_DEPS) 
	find test/ -name *out.txt -delete

$(TEST_TARGET) : $(TEST_OBJS) $(GCHS) $(TARGET) 
	$(CXXLINKER) $(LINKERFLAGS) -o $(TEST_TARGET) $(TEST_OBJS) $(TARGET)

test: $(TEST_TARGET)

gchs : $(GCHS)

info: 
	@echo sources:
	@echo $(SRCS)
	@echo deps:
	@echo $(DEPS)
	@echo headers:
	@echo $(HEADERS)
	@echo gchs:
	@echo $(GCHS)
	@echo test sources:
	@echo $(TEST_SRCS)

all:  info gchs $(TARGET) test
