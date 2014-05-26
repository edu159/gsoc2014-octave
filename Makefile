# Path definitions
ROOT_DIR=$(PWD)
BUILD_DIR=build
SRC_DIR=src


# Compiler flags
CPPFLAGS= -g

# Files to compile
OCT_FILES=ilu0.oct ilutp.oct

.SUFFIXES: .oct
.PHONY: benchmark

all : $(OCT_FILES)

%.oct : $(BUILD_DIR)/%.o
	mkoctfile $(CPPFLAGS) $<

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cc
	mkoctfile -c $(CPPFLAGS) $< -o $@

clean :
	rm *.oct

benchmark :
