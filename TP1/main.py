import sys
import numpy

# Função para ler a matriz blosum62 de um txt e retornar um dicionário para comparações
def get_blosum_dict(blosum62):
    lines = blosum62.readlines()
    blosum62.close()
    dic_blosum = {}
    x = lines[0]
    x = x.split()
    
    i = 1
    # Percorre o arquivo retornando as equivalências dos aminoacidos
    for i in range (1, (len(lines))):
        row = lines[i]
        row = row.split()
        j = 1
        for c in row[1:25]:
            dic_blosum[x[i-1],x[j-1]] = c 
            j+=1

    return(dic_blosum) 

# Função que cria a matriz de pontuações
def preenche_matriz(seq, blosum_dict):

	# Vê o tamanho das sequências para criar a matriz
	n = len(seq[0])+1	# adicionando v
	m = len(seq[1])+1	# adicionando w
	print(m, seq[0])
	matriz = numpy.zeros(shape=(n,m))

	# Faz as operações 
	penalidade = 0 


	# Preenche todas as posições da matriz com sua pontuação
	for i in range (1, n):
		for j in range (1, m):
			# Valor da BLOSUM62 para o casamento dos aminoácidos atuais da sequência 1 e 2
			blosum_score = (seq[0][i-1],seq[1][j-1])

			#print("Blosum:", blosum_score)
			#print("horizontal", matriz[i-1][j], int(blosum_dict.get(blosum_score)))
			#print("vertical", matriz[i][j-1], int(blosum_dict.get(blosum_score)))
			#print("diagonal", matriz[i-1][j-1], int(blosum_dict.get(blosum_score)))
			#print("")
			# Retorna o valor da posição horizontal anterior, o valor do casamento da BLOSUM62 para a posição e a penalidade usada para a execução
			horizontal = matriz[i-1][j]# + int(blosum_dict.get(blosum_score)) + penalidade
			# Retorna o valor da posição vertical anterior, o valor do casamento da BLOSUM62 para a posição e a penalidade usada para a execução
			vertical = matriz[i][j-1] #+ int(blosum_dict.get(blosum_score)) + penalidade

			# Se é um match na diagonal
			if(seq[0][i-1]==seq[1][j-1]):
				# Retorna o valor da posição diagonal anterior e o valor do casamento da BLOSUM62 para a posição
				diagonal = matriz[i-1][j-1] + int(blosum_dict.get(blosum_score))
			# Se mismatch na diagonal
			else:
				# Retorna o valor da posição diagonal anterior, o valor do casamento da BLOSUM62 para a posição e a penalidade usada para a execução
				diagonal = matriz[i-1][j-1] + int(blosum_dict.get(blosum_score)) + penalidade
			#print(horizontal, vertical, diagonal)

			# Pega o valor máximo dentre os três calculados para preencher a posição na matriz com as prioridades propostas em aula
			if (seq[0][i-1]==seq[1][j-1]):
				matriz[i][j]=diagonal
			elif (diagonal>=vertical and diagonal>=horizontal):
				matriz[i][j]=diagonal
			elif (horizontal>=vertical and horizontal>=horizontal):
				matriz[i][j]=horizontal
			else:
				matriz[i][j]=vertical
			#matriz[i][j]=max(horizontal, vertical, diagonal)


	return matriz

# Função que calcula o algoritmo Needleman Wunsch
def needleman_wunsch ():
	return

# Pega o nome do arquivo direto dos parâmetros de execução
file = sys.argv[1]

# Abre o arquivo para ler as sequências
f = open(file, 'r')

seq = []
seq = (f.read().splitlines() )

f.close()


blosum62 = open("BLOSUM62.txt", 'r')
blosum_dict = get_blosum_dict(blosum62)

#print(blosum_dict)
#print(len(blosum_dict))

matriz = preenche_matriz(seq, blosum_dict)
print(matriz)

