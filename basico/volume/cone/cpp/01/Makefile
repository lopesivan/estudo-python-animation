# Makefile Atualizado
V = 3.12.6
PY_PATH = /home/ivan/.pyenv/versions/$(V)
PY_SHORT_V = 3.12

# O caminho que você encontrou com o comando find (até a pasta anterior a /numpy)
NUMPY_INC = $(PY_PATH)/envs/jupyter/lib/python3.12/site-packages/numpy/_core/include

# Versão curta para diretórios internos
PY_SHORT_V = 3.12

CXXFLAGS = -std=c++17 -O3 \
           -I$(PY_PATH)/include/python$(PY_SHORT_V) \
           -I$(NUMPY_INC) \
           -I./matplotlib-cpp

# Flags do Linker
# Adicione -Wl,-rpath para que o binário saiba onde encontrar a .so correta em runtime
LDLIBS = -L$(PY_PATH)/lib -Wl,-rpath,$(PY_PATH)/lib -lpython$(PY_SHORT_V) -lcrypt -lpthread -ldl -lutil -lm# Makefile: Versão configurada para Python 3.12.6 com suporte a NumPy

all: volume_cone

volume_cone: volume_cone.o
	$(CXX) $(CXXFLAGS) volume_cone.o -o volume_cone $(LDLIBS)

clean:
	/bin/rm -f volume_cone volume_cone.o

# END OF FILE
