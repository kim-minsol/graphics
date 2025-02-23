Introduction to Computer Graphics (Spring 2023) - Assignment 6
=====
### Course/Assignment Skeleton Writers
Instructor: Minhyuk Sung (mhsung@kaist.ac.kr)

Last Update: May. 8, 2023 by Eunji Hong (hong969@kaist.ac.kr)

### Contents Writer
Minsol Kim (rlaalsthf02@naver.com)

# Keyframe List
이번 포스트는 Keyframe을 interpolation하여 animation을 제작하는 과정을 담았습니다. 

Keyframe animation을 구현하기 전에, keyframe을 저장하는 각각의 hotkey를 구현하였습니다.
- ‘u’: Update the current keyframe.
- ‘>’: Advance to the next keyframe.  
- ‘<’: Go back to the previous keyframe. 
- 'd’: Delete the current keyframe.  
- ‘n’: Create a new keyframe.  
- ' ': View the current keyframe.
- ‘i’: Import (read) a keyframe list file.  
- 'w’: Export (write) a keyframe list file.


Frame은 scene graph에 존재하는 transformation nodes의 RigTForms 집합으로 구성됩니다. 따라서, frame의 데이터 구조를 RigTForm 행렬의 vector로 표현할 수 있습니다.
'n'은 각각 현재 frame을 keyframe에 추가하는 기능을, 'u'는 현재 frame을 current key frame에 복사하는 기능을 맡기 때문에 RbtNode의 벡터인 Scene graph를 RigTForm 벡터 keyframe에 저장해야 합니다. 다른 모든 기능들은 반대로, keyframe을 Scene Graph에 전달하여 window로 확인할 수 있습니다.

---

# Linear Interpolation
가장 간단하게 두 프레임을 interpolation하는 방법은 선형적으로 구현하는 것입니다. 각각의 노드는 Translation, Rotation(Quaternion) vector를 갖고 있으므로, 두 벡터를 interpolation합니다.

## LERP
$$
\text{lerp}(\mathbf{c_0},\mathbf{c_1},\alpha)= (1-\alpha)\mathbf{c_0}+\alpha\mathbf{c_1}
$$   

Translation vector를 보간하는 것은 간단하게 두 벡터에 linear interpolation(LERP)을 적용하여 구현할 수 있습니다. 이때 $\alpha \in (0,1)$는 영상의 time에 따른 보간 상수가 되며, frame 간의 지정된 시간 ms를 영상 시간 ms에 나누어 계산됩니다.

## SLERP
![](https://i.imgur.com/mQEvAK2.png)
LERP는 Line을 따라 보간을 구현하므로 Translation에 적용이 가능합니다. 하지만, Rotation은 sphere에서의 보간이 구현되어야 합니다. 따라서 Spherical Linear Interpolation(SLERP)을 통해 회전 벡터를 보간합니다. Quaternion에 관한 설명은 포스트 아래에 따로 두었습니다.

## Quaternions Operations Overview
- scalar addition -> quaternion multiplication
- scalar negation -> quaternion inversion
- scalar multiplication -> quaternion power

## Power-based method
$$
\text{slerp}(\mathbf{q_0},\mathbf{q_{1}},\alpha) = (\text{cn}(\mathbf{q_1}\mathbf{q_{0}}^{-1}))^{\alpha}\mathbf{q_0}
$$   


- $\text{cn}$: conditional negate (첫 번째 좌표가 음수이면 negation 적용합니다.)   


이때 quaternions q와 -q는 같은 위치를 의미합니다. 따라서 더 짧은 보간을 위해 conditional negation을 적용하였습니다. 

## Quaternion method
$$
\frac{\sin((1-\alpha)\Omega)}{\sin(\Omega)} \overrightarrow{v}_{0} + \frac{\sin(\alpha\Omega)}{\sin(\Omega)} \overrightarrow{v}_{1}
$$   


n차원에 대해서 trigonometric argument(삼각함수 요소)들을 표현하면 SLERP를 다음과 같이 표현할 수 있습니다. 이때 $\Omega$는 벡터 사이의 각도를 의미합니다.

---

# Smooth Interpolation
Computer graphics에서 animation을 만들기 위한 기본적인 방법으로는 keyframe animation이 있습니다. discrete time에서의 keyframe들은 각각 다양한 objects와 camera의 locations, orientations 정보를 포함합니다. keyframe 사이를 smooth interpolation을 통해서 continuous 시간에 대해 parameter values를 채워 넣어 animation을 제작하는 방법입니다.
![](https://i.imgur.com/tmQJOBZ.png)
discrete values를 부드럽게 보간하는 간단한 방법은 *splines*을 이용하는 것입니다. Figure와 같이, 저차원의 polynomial function을 통해 독립된 값들(blue dots)을 연결시키는 함수가 splins가 됩니다.

## Cubic Bezier Functions
cubic polynomial function $c(t)$을 $t\in [0.,1]$에 대해 Beizer 표현으로 나타내는 method를 살펴보겠습니다. Beizer 표현에서는 $c_{0}, d_{0}, e_{0},c_{1}$ 총 4개의 control values로 cubic funtion을 표현합니다.

$$
\begin{aligned}
f &= (1-t)c_{0} + td_{0} \\
g &= (1-t)d_{0} + te_{0} \\
h &= (1-t)e_{0} + tc_{1} \\
m &= (1-t)f + tg \\
n &= (1-t)g+ th \\
c(t) &= (1-t)m + tn \\
\end{aligned}
$$   


해당 수식을 풀게 되면, 아래와 같이 간단하게 표현할 수 있습니다.   


$$
c(t)=c_{0}(1-t)^{3}+3d_{0}t(1-t)^{2}+3e_{0}t^{2}(1-t)+c_{1}t^{3}
$$   


cubic function에 0, 1의 값을 넣으면 $c(t)$와 control polygon이 각각 $c_{0}, c_{1}$에서 **접하는** 것을 확인할 수 있습니다. 모든 control value가 1인 경우 $c(t)=1$이 되므로, 단순 상수를 더하는 것이 됩니다.
또한, 해당 function에 미분을 가했을 때, $c'(0)=3(d_{0}-c_{0})$과 $c'(1)=3(c_{1}-e_{0})$의 값을 갖는 것을 확인할 수 있습니다. 따라서, $c'(i)=3(d_{i}-c_{i})=3(c_{i}-e_{i-1})$의 성질을 갖는 것도 확인할 수 있습니다.

### Translation
cubic function을 통해 $c_i$와 $c_{i+1}$ 두 값을 interpolation할 때의 다른 두 control points $d_{i},e_{i}$를 각각 구해야 합니다.

## Catmull-Rom Splines   
   
$$
\begin{aligned}
\frac{1}{2}(c_{i+1}-c_{i-1})=3(d_{i}-c_{i}) \rightarrow d_{i} &= \frac{1}{6}(c_{i+1}-c_{i-1}) + c_{i} \\
\frac{1}{2}(c_{i+2}-c_{i}) = 3(c_{i+1}-e_{i})
\rightarrow e_{i} &= \frac{-1}{6}(c_{i+2}-c_{i}) + c_{i+1} 
\end{aligned}
$$   

$d_{i},e_{i}$를 얻기 위해서, 이미 알고 있는 input data로부터의 $c_{i}$ 값들을 활용합니다. 더 정확히는 $c'(i)=3(d_{i}-c_{i})=3(c_{i}-e_{i-1})$의 성질을 이용하여 $d_{i}, e_{i}$의 값을 구하게 됩니다.
![](https://i.imgur.com/B1Pqcg8.png)   

$c'(t)=\frac{1}{2}(c_{i+1}-c_{i-1})$의 기울기와 $d_{i}$에 해당하는 점에서의 기울기 $c'(i)$가 같다는 점을,  $c'(t)=\frac{1}{2}(c_{i+2}-c_{i})$와 $e_i$에 해당하는 점에서의 기울기 $c'(i)$가 같다는 점을 이용하여 $d_i,e_i$를 구할 수 있습니다. 따라서, 각각 backward/forward로 $c_{-1}, c_{n+1}$의 값을 추가로 필요하게 됩니다.  위 Figure에서 화살표가 가르키는 부분이 같은 기울기를 갖는 부분입니다.

## Quaternion Splining   

$$
\begin{aligned}
d_{i}&=((c_{i+1}c^{-1}_{i-1})^{\frac{1}{6}})c_{i} 

\\

e_{i}&=((c_{i+2}c_{i}^{-1})^{\frac{-1}{6}})c_{i+1}
\end{aligned}
$$   

orientation에 대해 interpolation을 적용할 때, $x,y,z$에 대해 위의 방식을 직접적으로 적용하는 것은 불가능합니다. 따라서 간접적인 quaternion 연산으로 회전 보간을 구현할 수 있습니다. quaternion 연산으로 바꾸면 numerical 연산이 다음과 같이 바뀝니다.
- scalar add -> quaternion multiplication
- scalar negation -> quaternion inversion
- scalar multiplication -> quaternion power


따라서, 위의 고려사항을 Cubic Beizer function을 Quaternion에 적용하면, $d_{i},e_{i}$를 위와 같이 표현할 수 있습니다. 순차적인 alpha-blending을 적용하는 단계에서는 *lerp*를 *slerp*로 변경하여 구현할 수 있습니다.

---
# Quaternions   

$$
R = \begin{bmatrix} r  & 0\\ 0 & 1\end{bmatrix}
$$   

rotation matrix는 위와 같이 표현할 수 있으며, $[x,y,z]$에 대해 $\theta$ 각도로 회전시키는 행렬을 나타냅니다.   

$$
R_{\alpha}=(R_{1}R^{-1}_{0})^{\alpha}R_{0}
$$   


이를 통해 $\alpha$에 대한 보간을 진행하면 다음과 같은 matrix 보간이 가능합니다. $\alpha=0$일 때는 rotation matrix $R_{0}$, $\alpha=1$일 때는 rotation matrix $R_1$을 얻을 수 있습니다. 

이때의 문제점은 $R_{1}R^{-1}_{0}$과 같은 transition rotation matrix를 행렬의 axis/angle form에 맞게 고려해야 한다는 점입니다. quaternion을 사용하면 axis, angle에 대한 표현을 항상 갖기 때문에 추가적인 구현이 필요하지 않으면서도 rotation matrix의 역할을 수행할 수 있습니다. 또한, Rotation matrix에 비해 메모리가 적고, 연산도 효율적입니다.

## Quaternion representation   

$$
\begin{bmatrix} \cos(\frac{\theta}{2})\\\sin(\frac{\theta}{2})\hat{\mathbf{k}} \end{bmatrix}
$$   

Quaternion은 4개의 실수로 구성된 벡터입니다. 또한, axis $\hat{\mathbf{k}}$에 대해 $\theta$ 각도의 회전을 표현하기 위한 단위 벡터입니다. (norm == 1) 

## Operations
### Multiplication by a scalar   

$$
\alpha\begin{bmatrix} \cos(\frac{\theta}{2})\\\sin(\frac{\theta}{2})\hat{\mathbf{k}} \end{bmatrix} = \begin{bmatrix} \alpha\cos(\frac{\theta}{2}) \\ \alpha\sin(\frac{\theta}{2})\hat{\mathbf{k}} \end{bmatrix}
$$   


### Multiplication between two quaternions   

$$\begin{bmatrix} w_{1} \\ \hat{\mathbf{c}}_{1} \end{bmatrix}\begin{bmatrix} w_{2} \\ \hat{\mathbf{c}}_{2} \end{bmatrix} = \begin{bmatrix} (w_{1}w_{2}-\hat{\mathbf{c}}_{1}\cdot\hat{\mathbf{c}}_{2}) \\ w_{1}\hat{\mathbf{c}}_{2} + w_{2} \hat{\mathbf{c}}_{1} + \hat{\mathbf{c}}_{1}\times \hat{\mathbf{c}}_{2} \end{bmatrix}$$   


$\begin{bmatrix} w_{1} \\ \hat{\mathbf{c}}_{1} \end{bmatrix}$가 rotation matrix $R_1$을 표현한다고 가정하였을 때, 두 quaternion의 곱은 $R_{1}R_{2}$가 되는 특징으로 나타낼 수 있습니다.
![Gyiv9Us.png](https://i.imgur.com/Gyiv9Us.png)   

Quaternion의 vector 부분에 해당하는 i, j, k에 대한 곱 테이블입니다. 
![](https://i.imgur.com/gHpVECf.png)   

곱 테이블을 이용해서 두 사원수의 곱을 naive하게 계산하면 다음과 같고, 이를 간단하게 표현하면 $\begin{bmatrix} (w_{1}w_{2}-\hat{\mathbf{c}}_{1}\cdot\hat{\mathbf{c}}_{2})\\w_{1}\hat{\mathbf{c}}_{2}+w_{2}\hat{\mathbf{c}}_{1}+\hat{\mathbf{c}}_{1}\times \hat{\mathbf{c}}_{2} \end{bmatrix}$가 됩니다.

### Multiplicative inverse   

$$
\begin{bmatrix} \cos(\frac{\theta}{2})\\\sin(\frac{\theta}{2})\hat{\mathbf{k}} \end{bmatrix}^{-1} = \begin{bmatrix} \cos(\frac{\theta}{2})\\ -\sin(\frac{\theta}{2})\hat{\mathbf{k}} \end{bmatrix}
$$   

inverse 연산은 단순히 같은 축에서 $-\theta$로 회전하는 것이 됩니다.

### Quaternion Angle
quaternion angle을 계산할 때, arccos, arcsin만으로는 정확한 각도를 계산할 수 없습니다. $[0,2\pi]$일 때 $p=\cos\phi=\cos(-\phi), q=\sin(\phi)=\sin(\pi-\phi)$이므로 역삼각함수 적용하였을 때, unique angle을 얻을 수 없기 때문입니다. 따라서, arctan로 Quaternion의 각도를 계산합니다.

### Power   

$$
\begin{bmatrix} \cos(\frac{\theta}{2})\\\sin(\frac{\theta}{2})\hat{\mathbf{k}} \end{bmatrix}^{\alpha} = \begin{bmatrix} \cos(\frac{\alpha\theta}{2})\\ \sin(\frac{\alpha\theta}{2})\hat{\mathbf{k}} \end{bmatrix}
$$   

단위 축 $\hat{\mathbf{k}}$를 사원수로부터 추출한 후, arctan로 $\theta$를 추출합니다. 이후 $\alpha$를 곱하는 것으로 quaternion의 급수 연산을 구현할 수 있습니다.


# Result
![](https://i.imgur.com/tpV1Nrp.gif)





# Reference
[1] Foundations of 3D Computer graphics, [http://www.3dgraphicsfoundations.com/](http://www.3dgraphicsfoundations.com/)   

[2] OpenGL, [https://www.khronos.org/opengl/](https://www.khronos.org/opengl/)   

[3] Quaternions, Wikipedia, https://ko.wikipedia.org/wiki/%EC%82%AC%EC%9B%90%EC%88%98   

[4] Quaternion multiplication, https://gamesmith.tistory.com/161   

