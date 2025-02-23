Introduction to Computer Graphics (Spring 2023) - Assignment 8
=====
### Course/Assignment Skeleton Writers
Instructor: Minhyuk Sung (mhsung@kaist.ac.kr)

Last Update: May. 8, 2023 by Eunji Hong (hong969@kaist.ac.kr)

### Contents Writer
Minsol Kim (rlaalsthf02@naver.com)



# Contents

## Triangle Soup
OpenGL 렌더링에서 주요하게 사용되는 표현은 triangle soup입니다. 3개의 정점으로 표현되는 triangles의 집합으로 geometry를 표현합니다. 이때 geometry를 표현하기 위해, mesh를 사용합니다. (i.e. 어떤 삼각형들이 한 edge에서 만나는지 등)


## Mesh
![](https://i.imgur.com/jsivZBh.png)
mesh는 vertex, edge, face data로 표현되는 데이터 구조입니다. 위 Figure는 mesh data 구조에서 vertex, edge, face를 연결하는 과정의 시각화입니다.


## Types of coordinate function
- implicit function: $f(x,y,z)=0$
	- smooth surface를 표현하기 위해 surface points=0로 설정할 수 있습니다.
	- examples:
		- a plane: $x+2y-3z+1=0$
		- a sphere: $x^2+y^2+z^2-4=0$
- explicit function: $z=f(x,y)$
	- parametric representation: $(x(s,t),y(s,t),z(s,t))$

# Parametric Patches
## Tensor-product spline surface
**parametric patch** $(x(s,t),y(s,t),z(s,t))$를 통해 x,y,z에 대해 $(s,t)$ plane으로 표현할 수 있습니다. $(s,t)$ 평면에 대한 square, triangular 부분까지도 정의합니다. curves부터 surfaces까지 upgraded로 구현되며, m x n의 rectilinear control mesh로 폴리곤이 형성됩니다. s와 t 변수들에 대해 spline 정의를 적절하게 적용하여 parameteric patch를 완성합니다.



# Subdivision surfaces
![](https://i.imgur.com/FUwy06r.png)
patching 없이 single control mesh로 전체 geometry를 표현합니다. 하나의 메쉬에서 subdivision을 몇 번 적용하여 부드러운 표면을 근접하게 만들 수 있습니다. 이때 surface에 해당 과정을 몇 번만 적용합니다. 무한대로 적용하는 것도 가능하지만, 효율성을 위함입니다-control mesh로 부드러운 표면을 *근사*하는 것이 목적이기 때문입니다. 

## Catmull-Clark subdivision surface
watertight mesh $M^0$으로 시작하여 $M^1$으로 개선해가는 방식입니다.
- Connectivity of $M^{1}$:
	- $M^0$ 내 각 vertex $v$에 $M^1$에 위치하는 새로운 'vertex-vertex' $v_v$를 결합합니다.
	- $M^{0}$ 내 각 edge $e$에  $M^1$에 위치하는 새로운 'edge-vertex' $v_e$를 결합합니다.
	- $M^{0}$ 내 각 face $f$에  $M^1$에 위치하는 새로운 'face-vertex' $v_f$를 결합합니다.
위의 과정-위 Figure와 동일합니다.-으로 $M^{1}$을 설정하면, 모든 faces가 quads가 됩니다. 이때 $M^1$은 아래의 특징을 가집니다:
- valence의 모든 vertex가 four "ordinary"입니다.
- valence의 모든 vertex가 four "extraordinary"와 다릅니다.

### Terminology
- Ordinary vertices: Any vertex of valence four.
- Extraordinary vertices: All the rest.

### Geometry rules
![](https://i.imgur.com/WqYiEqk.png)

$M^i$에서의 face $f$가 vertices $v_{j}$에 둘러 쌓여 있을 때, 우리는 새로운 face-vertex $v_f$
를 아래의 식으로 생성합니다.

$$
v_{f}=\frac{1}{m_{f}}\sum\limits_{j}v_j
$$
$M^{i}$ vertices의 중점과 동일합니다. 위 식에서 vertices의 숫자 $m_{f}$은 4입니다.

### Edge rules
![](https://i.imgur.com/sckP09T.png)

다음은 $M^{i}$의 vertices $v_1, v_2$를 edge $e$로 잇는 과정입니다. 첫 번째로 faces $f_{1}$과 $f_2$를 분리한 후에, 새로운 $M^{i+1}$의 edge-vertex를 생성합니다.

$$
v_{e}=\frac{1}{4}(v_{1}+v_{2}+v_{f_{1}}+ v_{f_{2}})
$$
위의 수식처럼, edge를 구성하는 두 정점과 앞에서 새로 생성된 face-vertex 정점 간의 평균 값으로 구합니다.


### Vertex rules
![](https://i.imgur.com/IBGCYzp.png)

최종적으로 위 Figure에서의 v가 vertex가 되게 합니다. 앞의 두 과정은 새로운 정점을 생성했다면, vertex rules에서는 기존 정점들을 update합니다.
- $n_{v}$ vertices와 모두 연결되어 있습니다.
- $n_v$ faces $f_j$에 둘러 쌓여 있습니다.

$$
v_{v}=\frac{n_{v}-2}{n_{v}}v+ \frac{1}{n^{2}_{v}}\sum\limits_{j}v_{j}+ \frac{1}{n^{2}_{v}}\sum\limits_{j}v_{f_{j}}
$$
따라서, $M^{i+1}$의 새로운 vertex-vertex는 위와 같이 설정합니다. 

## Recursive subdivision
subdivision 과정을 $i \ge 1$에 대해 재귀적으로 반복합니다. 새로운 mesh가 생성되면 대략 4배의 vertex 수가 증가하는 것이 됩니다. subdivision을 반복할수록, 더 많은 mesh가 rectilinear하게 됩니다 --- extraordinary points는 이때 고정됩니다.


# Reference 
[1] https://mhsung.github.io/kaist-cs380-spring-2023/
[2] Steven J. Gortler, [Foundations of 3D Computer Graphics](https://mitpress.mit.edu/9780262017350/foundations-of-3d-computer-graphics/)
