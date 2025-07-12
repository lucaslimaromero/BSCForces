# Simulação de Feixes Acústicos

Este projeto é um conjunto de ferramentas desenvolvidas em Julia para simular e visualizar feixes de Bessel acústicos. O objetivo principal é construir a base para o cálculo de forças de radiação acústica em micropartículas.

## 💡 1. A Ideia Central 

O que acontece quando apontamos um feixe de som super focado para uma partícula minúscula? O som, sendo uma onda, carrega momento e pode literalmente empurrar, puxar ou até girar a partícula. Este projeto explora esse fenômeno.

Em vez de usar um feixe de som comum, que se espalha, usamos um **Feixe de Bessel**. A propriedade especial deste feixe é que ele mantém sua forma de "anéis concêntricos" por uma longa distância, quase como um "túnel de som".

Nosso objetivo é simular esse túnel de som e, futuramente, calcular a força que ele exerce.

## 🧩 2. Como Funciona? A Expansão em Ondas Parciais

Para descrever matematicamente um feixe tão complexo, usamos um método chamado **Expansão em Ondas Parciais**. A ideia é simples, e usa as duas vertentes Kantianas:

- **Síntese:** Quebramos o feixe de Bessel em uma soma de centenas de "ondinhas" esféricas muito simples.
- **Análise:** Para cada "ondinha" de ordem `(n, m)`, precisamos de um peso exato para que a soma final reconstrua o feixe de Bessel. Essa "quantidade" é um número complexo chamado **Coeficiente de Forma do Feixe (BSC)**, ou `A_n^m`.

O nosso código, portanto, faz duas coisas essenciais:
1.  **Calcula a "Receita":** Determina os BSCs para o feixe que queremos simular.
2.  **Monta o Feixe:** Soma todas as "ondinhas" para visualizar o campo de intensidade em qualquer ponto do espaço.

## 📜 3. Estrutura do Código

O projeto é organizado em torno de algumas funções chave:

#### `bsccalc(n, m, α)`
- **Responsabilidade:** Calcular o Coeficiente de Forma do Feixe (`A_n^m`) para um feixe **centrado na origem**.
- **Entradas:** Índices `n` e `m` da onda parcial, e o ângulo do áxicon `α` do feixe.
- **Saída:** O valor do BSC para aquela componente.

#### `partialwavexp(...)`
- **Responsabilidade:** Calcular a amplitude total da onda em um único ponto do espaço.
- **Como funciona:** Recebe as coordenadas do ponto `(r, θ, ϕ)` e as "receitas" (BSCs e Polinômios de Legendre) já pré-calculadas para ser extremamente rápido. Ela realiza o somatório das ondas parciais.

#### `makeplotpartial(...)` e `makeplotpartial_displaced(...)`
- **Responsabilidade:** Orquestrar a simulação e gerar os gráficos de intensidade.
- **Como funcionam:**
    1. Definem um grid de pixels para a imagem.
    2. Realizam o **pré-cálculo** dos BSCs e dos Polinômios de Legendre (a etapa mais importante para a performance).
    3. Para cada pixel do grid, chamam a `partialwavexp` para calcular a intensidade.
    4. Usam a biblioteca `Plots.jl` para montar e salvar a imagem final.
- A versão `_displaced` implementa o truque de **transladar o sistema de coordenadas** (`x_rel = xi - x₀`, etc.) para simular um feixe deslocado de forma eficiente, sem precisar de um somatório duplo.

## ▶️ 4. Como Usar 

1.  **Inicie o Julia com Threads:** Para máxima velocidade, use o paralelismo.
    ```bash
    julia -t auto
    ```
2.  **Execute o Script:** Dentro do REPL do Julia, inclua e execute o arquivo principal.
    ```julia
    include("seu_script.jl")
    ```
3.  **Analise os Resultados:** O script irá gerar as imagens `.png` na pasta do projeto. Modifique os parâmetros no final do arquivo (resolução da imagem, deslocamento, `n_max`) para explorar diferentes configurações.
