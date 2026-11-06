# Asignación presupuestaria con Gradientes Conjugados Precondicionado (GCP)

- **Estudiante:** Camila Sarahí Arias Vargas
- **Mención:** Desarrollo de Sistemas
- **Área:** Finanzas
- **Página web:** https://sarav121.github.io/GCP_pagina/asignacion.html

## Descripción

Simulador web que reparte un presupuesto entre ocho categorías de gasto (Personal, Infraestructura, Marketing digital, I+D, Capacitación, Seguridad informática, Soporte y operaciones, Legal y cumplimiento). El problema se formula como un sistema lineal **Ax = b**, con **A** simétrica y definida positiva, y se resuelve con:

1. un **método directo** (factorización de Cholesky), usado como solución numérica de referencia;
2. el **Método de Gradientes Conjugados Precondicionado (GCP)**, basado en Suñagua (2020);
3. el **Gradiente Conjugado simple (GC)**, solo como contraste.

El análisis se hace en tres escenarios: **Ideal**, **Bajo estrés** y **Mal condicionado**.

## Modelo financiero

Se busca la asignación **x** (miles de Bs) que minimiza el costo de desajuste

f(x) = ½ xᵀAx − bᵀx

Como A es simétrica y definida positiva, el mínimo es único y cumple Ax = b.

| Variable | Categoría |
|---|---|
| x₁ | Personal |
| x₂ | Infraestructura |
| x₃ | Marketing digital |
| x₄ | I+D |
| x₅ | Capacitación |
| x₆ | Seguridad informática |
| x₇ | Soporte y operaciones |
| x₈ | Legal y cumplimiento |

- **aᵢᵢ:** costo marginal de desajustar la categoría i.
- **aᵢⱼ (i ≠ j):** acoplamiento entre categorías.
- **bᵢ:** necesidad prioritaria de la categoría i.

Los datos son ilustrativos y están definidos para este ejercicio académico.

## Escenarios

| Escenario | Qué cambia | Precondicionador por defecto |
|---|---|---|
| Ideal | A con diagonal dominante y poco acoplamiento | Jacobi |
| Bajo estrés | Inflación: diagonal de A +10 % y aumento de b | Jacobi |
| Mal condicionado | A = DCD con una matriz diagonal D de escalas muy distintas: mal escalamiento diagonal | Jacobi |

## Funcionalidades

- Edición de los coeficientes de A, del vector b, del vector inicial x₀, de la tolerancia y del máximo de iteraciones, con recálculo automático.
- Selección del precondicionador: Jacobi, SSOR (ω = 1) o ninguno.
- Validación de A: simetría, menores principales (Sylvester) y autovalores. Si A no es válida, se muestra un aviso y no se calcula.
- Solución de referencia, GCP vs GC, errores por variable (absoluto, relativo y porcentual), residuo relativo global y κ(A), κ(M⁻¹A).
- Gráfico de convergencia del residuo relativo.
- Tabla comparativa de los tres escenarios.
- Módulo de pruebas QA con formato Entrada, Resultado esperado, Resultado obtenido y Conclusión.

## Estructura del repositorio

```
.
├── asignacion.html   # Estructura y contenido de la página
├── style.css   # Estilos (modo claro y oscuro)
├── script.js       # Lógica numérica, validación, interfaz y pruebas QA
└── README.md
```

## Arquitectura

Es una aplicación de una sola página, sin frameworks ni dependencias, salvo las fuentes de Google. Todo el cálculo se hace en el navegador.

**Capa numérica (`script.js`)**

| Función | Rol |
|---|---|
| `chol`, `fsub`, `bsub`, `csolve` | Factorización de Cholesky y sustituciones; solución directa |
| `det` | Determinantes de los menores principales (Sylvester) |
| `eig`, `kappa` | Autovalores por rotaciones de Jacobi; κ = λmax / λmin |
| `validate` | Comprueba datos numéricos, simetría y definición positiva |
| `precM` | Construye la matriz precondicionadora M (Jacobi, SSOR o identidad) |
| `pcg` | Algoritmo GCP; con M = I equivale al GC simple |
| `kPre` | κ(M⁻¹A) mediante la matriz simétrica L⁻¹AL⁻ᵀ, con M = LLᵀ |
| `analyze` | Ejecuta referencia, GCP, GC y condicionamiento de un escenario |

**Capa de interfaz (`script.js`)**

| Función | Rol |
|---|---|
| `build`, `load`, `read` | Crea las entradas, carga un escenario y lee los datos de pantalla |
| `render` | Recalcula y dibuja todas las secciones |
| `cmp`, `chart` | Tabla comparativa de escenarios y gráfico de convergencia |
| `qa` | Ejecuta y muestra las pruebas de software |

## Algoritmo GCP

Con x₀ (por defecto 0) y r₀ = b − Ax₀:

1. Resolver M z₀ = r₀ y hacer p₀ = z₀.
2. Para k = 0, 1, 2, …:
   - αₖ = (rₖ · zₖ) / (pₖ · A pₖ)
   - xₖ₊₁ = xₖ + αₖ pₖ
   - rₖ₊₁ = rₖ − αₖ A pₖ
   - Resolver M zₖ₊₁ = rₖ₊₁
   - βₖ = (rₖ₊₁ · zₖ₊₁) / (rₖ · zₖ)
   - pₖ₊₁ = zₖ₊₁ + βₖ pₖ

**Criterio de parada:** ‖rₖ‖₂ / ‖b‖₂ < 10⁻⁶, o alcanzar el máximo de iteraciones.

## Pruebas de software (QA)

| Prueba | Qué verifica |
|---|---|
| Carga de escenarios | Los tres escenarios son válidos y GCP converge |
| Matriz no válida | Se rechazan una matriz no simétrica y una no definida positiva, sin colapsar |
| Edición dinámica | Al cambiar un coeficiente, la solución se recalcula y GCP coincide con la referencia |
| Precondicionamiento | En el escenario mal condicionado, GCP necesita ≤ iteraciones que GC |
| Validación de datos actuales | Comprueba la matriz y el vector que están cargados en la interfaz |

Se ejecutan desde la sección **Módulo de pruebas (QA)** con el botón **Ejecutar pruebas**.

## Limitaciones

- La interfaz trabaja con n = 8 variables (ocho categorías de gasto).
- No está implementado el precondicionador de Cholesky incompleto (IC0); se ofrecen Jacobi y SSOR.
- En aritmética exacta, el GC termina en n iteraciones como máximo; con n = 8 la ventaja del precondicionador puede verse en el número de iteraciones y en κ(M⁻¹A).

## Referencia


Suñagua S., P. (2020). *Método de Gradientes Conjugados Precondicionado*. Revista Boliviana de Matemática, UMSA.
