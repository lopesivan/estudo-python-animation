# Explicação Matemática: Método das Cascas Cilíndricas (Cálculo de Volume)

Vou explicar com cuidado conceitual e matemático, passo a passo, justificando a fórmula do $dV$ usada no método das cascas no seu código.

### Observação importante (conceitual):
No seu exemplo do cone, o método correto é o método das **cascas cilíndricas**, não esféricas.
* “Cascas esféricas” são usadas quando o sólido é naturalmente descrito por esferas (ex.: uma esfera carregada).
* Aqui, como a variável é $\rho$ (raio no plano), estamos lidando com cascas cilíndricas.

---

### 1. Ideia geométrica do método das cascas cilíndricas
Em vez de decompor o cone em discos horizontais, nós o decompomos em cascas verticais:
* Cada casca é um cilindro fino e oco;
* **Raio:** $\rho$
* **Espessura:** $d\rho$
* **Altura:** depende de $\rho$

**Visualmente:** É como “descascar” o cone em tubos concêntricos, do eixo até a borda.

### 2. Volume elementar de uma casca cilíndrica
O volume de uma casca cilíndrica fina é:

$$\boxed{dV = (\text{área lateral}) \times (\text{espessura})}$$

A área lateral de um cilindro de raio $\rho$ e altura $h(\rho)$ é:

$$A_{\text{lateral}} = 2\pi\rho \cdot h(\rho)$$

Logo:

$$\boxed{dV = 2\pi\rho \cdot h(\rho) \, d\rho}$$

### 3. Determinação da altura da casca $h(\rho)$
O cone tem:
* **Raio da base:** $r$
* **Altura total:** $h$

A altura do cone decresce linearmente com o raio $\rho$. Isso vem de semelhança de triângulos:
* Quando $\rho = 0$, a altura é $h$.
* Quando $\rho = r$, a altura é $0$.

Relação geométrica:
$$\frac{h(\rho)}{h} = 1 - \frac{\rho}{r}$$

Portanto:
$$\boxed{h(\rho) = h\left(1 - \frac{\rho}{r}\right)}$$
ou, expandindo:
$$\boxed{h(\rho) = h - \frac{h}{r}\rho}$$

### 4. Substituindo na expressão de $dV$
Voltamos à fórmula geral da casca:
$$dV = 2\pi\rho \cdot h(\rho) \, d\rho$$

Substituindo $h(\rho)$:
$$dV = 2\pi\rho \left(h - \frac{h}{r}\rho\right) d\rho$$

E chegamos exatamente à fórmula que você usou:
$$\boxed{dV = 2\pi\rho \left(h - \frac{h}{r}\rho\right) d\rho}$$

### 5. Por que os limites são $0 \le \rho \le r$
Como estamos varrendo as cascas:
* Começamos no eixo do cone: $\rho = 0$
* Terminamos na borda da base: $\rho = r$

Logo, o volume total é:
$$V = \int_0^r 2\pi\rho\left(h-\frac{h}{r}\rho\right)d\rho$$

### 6. Resultado final (consistência física)
Ao resolver a integral, obtemos:
$$V = \frac{1}{3}\pi r^2 h$$

Isso confirma a coerência geométrica e a consistência com o volume conhecido do cone.

### 7. Resumo conceitual
1.  O método das cascas cilíndricas decompõe o sólido em tubos finos.
2.  Cada tubo contribui com $dV = (\text{perímetro}) \times (\text{altura}) \times (\text{espessura})$.
3.  No cone, a altura varia linearmente com o raio.
