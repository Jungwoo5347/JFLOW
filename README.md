# JFLOW

2차원 압축성 오일러 방정식을 푸는 유한체적(FVM) 솔버입니다. Fortran으로 작성했습니다.
곡선 정렬 격자(structured curvilinear grid)에서 셀 중심(cell-centered) 방식으로 계산하고,
Pointwise 등에서 만든 PLOT3D 격자(`.x`)를 그대로 읽습니다.

![RAE 2822 transonic airfoil](docs/images/rae2822_cp.png)

## 기능

| 구분 | 내용 |
|---|---|
| 지배 방정식 | 2D 압축성 오일러 방정식 (비점성, 이상기체) |
| 격자 | 단일 블록 곡선 정렬 격자. PLOT3D ASCII / 바이너리 (2D, 3D 의 k=1 면), 균일 직교 격자 자동 생성 |
| 리만 솔버 | Rusanov, Roe |
| 공간 정확도 | 1차, 또는 MUSCL κ 기법 + minmod 리미터 (특성변수 재구성, 면 값이 음수면 1차로 되돌림) |
| 시간 적분 | SSP-RK3 (3단계 TVD Runge–Kutta) |
| 정상 문제 | 국소 시간 간격 (셀마다 자기 Δt), 정규화 잔차 L2 norm 수렴 판정 |
| 경계조건 | 벽(미끄럼), 대칭면, 초음속 입구, 출구, 리만 불변량 원방 |
| 격자 연결 | C-격자 wake cut, O-격자 이음매를 자동으로 찾아 연결 |
| 병렬화 | OpenMP |

## 빌드

gfortran이 필요합니다 (Windows에서 테스트).

```
build.bat                   디버그 빌드 (배열 범위 검사, backtrace)
build.bat release           최적화 빌드 (-O3)
build.bat release omp       최적화 + OpenMP
```

빌드가 끝나면 `main.exe` 가 만들어집니다.

## 실행

케이스 폴더 하나에 입력 파일 두 개(`grid.in`, `flow.in`)를 넣고, 그 폴더를 인자로 줍니다.

```powershell
$env:OMP_NUM_THREADS = 8        # OpenMP 빌드일 때 스레드 수
.\main.exe cases\sod
```

결과는 같은 폴더에 저장됩니다. 파일 이름 앞부분은 `grid.in` 의 `$output` 값입니다.

| 파일 | 내용 |
|---|---|
| `<output>_result.dat` | 전체 해 (Tecplot ASCII, 노드 좌표 + 셀 중심 ρ, u, v, p, Mach) |
| `<output>_line.dat` | 셀 한 줄의 값 (`$line_row`), 1D 그래프용 |
| `<output>_residual.dat` | 매 반복의 정규화 잔차 (연속, x/y 운동량, 에너지) |

## 예제: Sod 충격관

격자 파일 없이 바로 돌릴 수 있는 예제입니다. 폴더(예: `cases/sod`)를 만들고 아래 두 파일을 넣으세요.

`grid.in`
```
$block-1 = cartesian 101 21 0.0 1.0 0.0 0.2     # 균일 격자: ni nj xmin xmax ymin ymax
$bc-1-1  = 220
$bc-1-2  = 220
$bc-1-3  = 220
$bc-1-4  = 220
$output  = sod
```

`flow.in`
```
$steady         = 0              # 비정상 문제
$t_end          = 0.2
$cfl            = 0.4
$print_interval = 10
$line_row       = 10             # sod_line.dat 로 뽑을 셀 줄
$gamma          = 1.4

$sigma          = 1              # MUSCL
$kappa          = 0.3333333333333333
$b              = 1
$flux           = roe

$state-left     = 1.0    0.0  0.0  1.0     # rho u v p
$state-right    = 0.125  0.0  0.0  0.1
$initial        = right
$region-1       = left   *  0.5   *  *     # x < 0.5 인 셀은 left 상태
```

```powershell
.\main.exe cases\sod
```

`cases/sod/sod_line.dat` 의 열은 x, ρ, u, v, p 입니다. t = 0.2 의 정확해와 비교하면 됩니다.

## 입력 파일

두 파일 모두 `$키 = 값` 형식이고, `#` 뒤는 주석입니다 (KFLOW `grid.in` 형식과 호환).
전체 키 설명은 [case_setup.f90](case_setup.f90) 맨 위에 있습니다.

### grid.in

| 키 | 설명 |
|---|---|
| `$block-1` | 격자 파일 경로 (케이스 폴더 기준), 또는 `cartesian ni nj xmin xmax ymin ymax` |
| `$form-1` | 0: 바이너리 (unformatted), 1: ASCII |
| `$bc-1-f = 코드 [n1 n2] [상태]` | 면 f (1 imin, 2 imax, 3 jmin, 4 jmax) 의 노드 n1~n2 에 경계조건. 범위를 생략하면 면 전체, 뒤에 쓴 줄이 앞의 줄을 덮어씀 |
| `$output` | 결과 파일 이름 앞부분 |

경계조건 코드

| 코드 | 경계조건 |
|---|---|
| 110 | 벽 (비점성 미끄럼) |
| 200 | 원방 (리만 불변량, 아음속/초음속 유입·유출 자동 처리) |
| 210 | 초음속 입구 (상태 고정) |
| 220 | 출구 (안쪽 값 외삽) |
| 300 | 대칭면 |

경계조건을 주지 않은 면은 양 끝 노드가 겹치는 다른 면을 찾아 자동으로 연결합니다 (C-격자 wake cut, O-격자 이음매).
짝을 찾지 못한 면이 남으면 위치를 알려주고 멈춥니다.

### flow.in

| 키 | 설명 |
|---|---|
| `$steady` | 1: 정상 (수렴할 때까지), 0: 비정상 (`$t_end` 까지) |
| `$max_iter`, `$res_tol` | 정상 문제의 최대 반복, 수렴 기준 (`0` 이면 판정 없이 `max_iter` 까지) |
| `$local_dt` | 1: 국소 시간 간격 (정상 문제만), 0: 전체 최소 Δt |
| `$cfl`, `$gamma` | CFL 수, 비열비 |
| `$sigma`, `$kappa`, `$b` | MUSCL on/off, κ, minmod 압축 계수 (1 ≤ b ≤ (3−κ)/(1−κ)) |
| `$flux` | `rusanov` 또는 `roe` |
| `$state-이름 = ρ u v p` | 이름 붙인 유동 상태 (경계조건, 초기조건에서 사용) |
| `$initial`, `$region-k` | 초기 상태, 상자 영역만 다른 상태로 (`*` 는 무한대) |

무차원화는 ρ∞ = 1, p∞ = 1/γ 로 잡으면 음속이 1이 되어 속도가 곧 마하수가 됩니다.
받음각 α 는 `u = M cos α`, `v = M sin α` 로 줍니다.

## 결과

**Sod 충격관** (100셀): 성분별(원시변수) 재구성은 접촉면 뒤에 내부에너지 오버슈트가 생기고,
특성변수 재구성은 이를 없앱니다.

![Sod shock tube](docs/images/sod_contact_shock.png)

**다이아몬드 에어포일** (M 2.0, α 2°): 충격파–팽창파 이론과 비교. 면별 압력 오차 1.5% 이내.

![Diamond airfoil](docs/images/diamond_airfoil.png)

**RAE 2822 천음속 에어포일** (M 0.729, α 2.52°): Airbus pETW 풍동 시험 Cp
([NWTF 데이터셋](https://www.nwtf.ac.uk/dataset/2687/))과 비교.
비점성 계산이라 경계층 효과가 없어서 충격파가 조금 뒤에, 더 강하게 서고 양력이 크게 나옵니다.
(맨 위 그림: 왼쪽 표면 Cp, 오른쪽 마하수 분포. 위 면에 국소 초음속 영역과 충격파가 보입니다)

**과소팽창 제트** (노즐 M 2.0, 압력비 1.5, 주변류 M 1.5): 충격 다이아몬드 구조.

![Underexpanded jet](docs/images/underexpanded_jet.png)

**국소 시간 간격의 수렴 가속** (원기둥 M 3, SC1095 에어포일):

![Local time step](docs/images/local_time_step.png)

## 소스 구조

| 파일 | 내용 |
|---|---|
| [euler.f90](euler.f90) | 상수, 보존변수 ↔ 원시변수 변환, 음속 |
| [case_setup.f90](case_setup.f90) | `grid.in`, `flow.in` 읽기 |
| [grid.f90](grid.f90) | 격자 생성 / PLOT3D 읽기, 셀 면적과 면 법선 계산 |
| [flow.f90](flow.f90) | 유동 배열, 초기조건, 경계조건 (고스트 셀 2층), 자동 면 연결 |
| [solver.f90](solver.f90) | MUSCL 재구성, 리만 솔버, 잔차, 시간 간격 |
| [monitor.f90](monitor.f90) | 잔차 기록과 수렴 판정, 발산 검사 |
| [output.f90](output.f90) | Tecplot 출력 |
| [main.f90](main.f90) | 입력 → 격자 → 시간 루프 (SSP-RK3) → 출력 |
