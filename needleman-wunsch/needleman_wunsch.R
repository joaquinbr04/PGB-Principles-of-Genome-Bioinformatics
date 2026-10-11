## STEP 1: CREAR LA MATRIZ DE PUNTUACIONES

blossum50 <- read.table(file = "blossum50.txt")

seq1 <- "BEAVTY"
seq2 <- "BEAST"

# Convertimos cada secuencia en un vector de aminoácidos
aa1 <- strsplit(seq1, "")[[1]]
aa2 <- strsplit(seq2, "")[[1]]

seq1length <- length(aa1)
seq2length <- length(aa2)

# Las filas representan seq2 y las columnas representan seq1
m <- matrix(
  NA,
  nrow = seq2length + 1,
  ncol = seq1length + 1
)

rownames(m) <- c("-", aa2)
colnames(m) <- c("-", aa1)

# Primera fila y primera columna
m[, 1] <- 0
m[1, ] <- 0

# Rellenamos las puntuaciones BLOSUM50
for (i in 2:nrow(m)) {
  for (j in 2:ncol(m)) {
    
    m[i, j] <- blossum50[
      aa2[i - 1],
      aa1[j - 1]
    ]
  }
}


## STEP 2: MATRIZ DE PUNTUACIONES ACUMULADAS

gap <- -1

AS <- matrix(
  NA,
  nrow = nrow(m),
  ncol = ncol(m)
)

rownames(AS) <- rownames(m)
colnames(AS) <- colnames(m)

# Matriz de direcciones
direction <- matrix(
  "",
  nrow = nrow(m),
  ncol = ncol(m)
)

rownames(direction) <- rownames(m)
colnames(direction) <- colnames(m)

# Esquina inicial
AS[1, 1] <- 0
direction[1, 1] <- "Done"

# Inicializamos la primera fila
for (j in 2:ncol(AS)) {
  AS[1, j] <- AS[1, j - 1] + gap
  direction[1, j] <- "L"
}

# Inicializamos la primera columna
for (i in 2:nrow(AS)) {
  AS[i, 1] <- AS[i - 1, 1] + gap
  direction[i, 1] <- "U"
}

# Rellenamos las celdas interiores
for (i in 2:nrow(AS)) {
  for (j in 2:ncol(AS)) {
    
    diagonal <- AS[i - 1, j - 1] + m[i, j]
    up       <- AS[i - 1, j] + gap
    left     <- AS[i, j - 1] + gap
    
    best <- max(diagonal, up, left)
    
    AS[i, j] <- best
    
    moves <- character(0)
    
    if (diagonal == best) {
      moves <- c(moves, "D")
    }
    
    if (up == best) {
      moves <- c(moves, "U")
    }
    
    if (left == best) {
      moves <- c(moves, "L")
    }
    
    direction[i, j] <- paste(moves, collapse = "/")
  }
}


## STEP 3: RECONSTRUIR EL ALINEAMIENTO

# Empezamos en la esquina inferior derecha
i <- nrow(direction)
j <- ncol(direction)

aligned1 <- character(0)
aligned2 <- character(0)

# Continuamos hasta que i y j sean ambos iguales a 1
while (i > 1 || j > 1) {
  
  move <- direction[i, j]
  
  # Movimiento diagonal
  if (grepl("D", move)) {
    
    aligned1 <- c(aligned1, aa1[j - 1])
    aligned2 <- c(aligned2, aa2[i - 1])
    
    i <- i - 1
    j <- j - 1
    
    # Movimiento hacia arriba
  } else if (grepl("U", move)) {
    
    aligned1 <- c(aligned1, "-")
    aligned2 <- c(aligned2, aa2[i - 1])
    
    i <- i - 1
    
    # Movimiento hacia la izquierda
  } else if (grepl("L", move)) {
    
    aligned1 <- c(aligned1, aa1[j - 1])
    aligned2 <- c(aligned2, "-")
    
    j <- j - 1
  }
}

# Invertimos porque hemos recorrido las matrices al revés
aligned1 <- rev(aligned1)
aligned2 <- rev(aligned2)

# Convertimos los vectores en texto
alignment1 <- paste(aligned1, collapse = "")
alignment2 <- paste(aligned2, collapse = "")

# Mostramos el resultado
cat(alignment1, "\n")
cat(alignment2, "\n")