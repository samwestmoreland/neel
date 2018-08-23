#
# Makefile for neighbour list code
#
# Makefile inspired by Job Vranish at
# https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
#

TARGET_EXEC ?= neel

GCC = g++-8

BUILD_DIR ?= ./build
SRC_DIRS ?= src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLGS) -MMD -MP -std=c++11

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(GCC) $(OBJS) -o $@ $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(GCC) $(CPPFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p
