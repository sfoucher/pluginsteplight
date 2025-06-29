# STL_Grid3D Processing Workflows

This document contains Mermaid flowcharts explaining the various processing workflows for STL_Grid3D objects.

## 1. Main Grid Creation Workflow

```mermaid
flowchart TD
    A[Point Cloud Input] --> B[Calculate Bounding Box]
    B --> C[Set Grid Resolution]
    C --> D[Create Empty Grid]
    D --> E[Initialize Ray Length Grid]
    E --> F[Multi-threaded Processing]
    
    F --> G[Thread 1: Process Points 0-N]
    F --> H[Thread 2: Process Points N-2N]
    F --> I[Thread N: Process Points N*M-End]
    
    G --> J[Create STL_Grid3D<int>]
    H --> K[Create STL_Grid3D<int>]
    I --> L[Create STL_Grid3D<int>]
    
    J --> M[Ray Tracing Algorithm]
    K --> M
    L --> M
    
    M --> N[STL_Grid3DBeamVisitor]
    N --> O[STL_Grid3DWooTraversalAlgorithm]
    O --> P[Accumulate Values in Grid]
    
    P --> Q[Merge Thread Results]
    Q --> R[Set Point Cloud Pointers]
    R --> S[Set Bounding Box]
    S --> T[Attach Ray Length Grid]
    T --> U[Compute Min/Max Values]
    U --> V[Output STL_Grid3D]
```

## 2. Ray Tracing Algorithm Workflow

```mermaid
flowchart TD
    A[Input: Point + Normal] --> B[Create CT_Beam]
    B --> C[Calculate End Points]
    C --> D[Initialize Woo Traversal]
    D --> E[Set Grid Boundaries]
    E --> F[Calculate Initial tMax Values]
    F --> G[Calculate DeltaT Values]
    
    G --> H{Keep First Cell?}
    H -->|Yes| I[Visit Current Cell]
    H -->|No| J[Skip Current Cell]
    
    I --> K[Add Value to Grid]
    J --> L[Calculate Next Step Axis]
    K --> L
    
    L --> M[Move to Next Voxel]
    M --> N{Within Grid Bounds?}
    N -->|No| O[End Traversal]
    N -->|Yes| P[Get Cell Index]
    P --> Q[Visit Cell with Visitors]
    Q --> R[Update tMax Values]
    R --> L
    
    O --> S[Return Accumulated Values]
```

## 3. Grid Filtering Workflows

### 3.1 Neighbor-Based Filtering

```mermaid
flowchart TD
    A[Input STL_Grid3D] --> B[Create Filtered Grid Copy]
    B --> C[Initialize Loop: Z=0 to Zmax]
    C --> D[Initialize Loop: X=0 to Xmax]
    D --> E[Initialize Loop: Y=0 to Ymax]
    
    E --> F[Get Current Cell Value]
    F --> G{Value == 0?}
    G -->|Yes| H[Skip to Next Cell]
    G -->|No| I[Set is_maxima = true]
    
    I --> J[Check Neighbors Loop]
    J --> K[Calculate Neighbor Position]
    K --> L{Neighbor in Grid Bounds?}
    L -->|No| M[Skip Neighbor]
    L -->|Yes| N[Get Neighbor Value]
    
    N --> O{Current < Neighbor?}
    O -->|Yes| P[Set is_maxima = false]
    O -->|No| Q[Continue Checking]
    
    P --> R[Break Neighbor Loop]
    Q --> J
    M --> J
    
    R --> S{is_maxima == false?}
    S -->|Yes| T[Set Cell to 0]
    S -->|No| U[Keep Cell Value]
    
    T --> H
    U --> H
    H --> V{More Cells?}
    V -->|Yes| E
    V -->|No| W[Return Filtered Grid]
```

### 3.2 Threshold-Based Filtering

```mermaid
flowchart TD
    A[Input STL_Grid3D] --> B[Create Filtered Grid Copy]
    B --> C[Set Threshold Value]
    C --> D[Initialize Iterator]
    D --> E{More Pixels?}
    E -->|No| F[Return Filtered Grid]
    E -->|Yes| G[Get Current Pixel Value]
    G --> H{Value < Threshold?}
    H -->|Yes| I[Set Pixel to 0]
    H -->|No| J[Keep Pixel Value]
    I --> K[Next Pixel]
    J --> K
    K --> E
```

### 3.3 Fast Filter with Ray Tracing

```mermaid
flowchart TD
    A[Input STL_Grid3D] --> B[Create Working Grid Copy]
    B --> C[Create Filtered Grid Copy]
    C --> D[Create Filter Visitor]
    D --> E[Create Set Value Visitor]
    E --> F[Initialize Traversal Algorithms]
    
    F --> G[Loop Through Point Cloud]
    G --> H[Get Point + Normal]
    H --> I[Create Beam 1: Point + Normal]
    I --> J[Create Beam 2: Point - Normal]
    
    J --> K[Set End Points]
    K --> L[Reset Filter Visitor]
    L --> M[Traverse Beam 1]
    M --> N[Get Votes for Beam 1]
    
    N --> O[Reset Filter Visitor]
    O --> P[Traverse Beam 2]
    P --> Q[Get Votes for Beam 2]
    
    Q --> R[Calculate Ratio: Max/Min]
    R --> S{Ratio < Threshold?}
    S -->|Yes| T[Filter Both Directions]
    S -->|No| U{Beam 1 > Beam 2?}
    
    U -->|Yes| V[Filter Beam 2 Only]
    U -->|No| W[Filter Beam 1 Only]
    
    T --> X[Set Zero in Both Directions]
    V --> Y[Set Zero in Beam 2 Direction]
    W --> Z[Set Zero in Beam 1 Direction]
    
    X --> AA{More Points?}
    Y --> AA
    Z --> AA
    
    AA -->|Yes| G
    AA -->|No| BB[Compute Min/Max]
    BB --> CC[Return Filtered Grid]
```

## 4. Local Maxima Detection Workflow

```mermaid
flowchart TD
    A[Input STL_Grid3D] --> B[Set Neighborhood Size]
    B --> C[Clear Output Vector]
    C --> D[Initialize Loop: X=0 to Xmax]
    D --> E[Initialize Loop: Y=0 to Ymax]
    E --> F[Initialize Loop: Z=0 to Zmax]
    
    F --> G[Get Current Cell Value]
    G --> H{Value > 0?}
    H -->|No| I[Skip to Next Cell]
    H -->|Yes| J[Create Pixel Position]
    
    J --> K{Within BBox Constraints?}
    K -->|No| I
    K -->|Yes| L[Check if Local Maxima]
    
    L --> M[Get Neighborhood Bounds]
    M --> N[Check All Neighbors]
    N --> O{Current > All Neighbors?}
    O -->|Yes| P[Add to Local Maxima List]
    O -->|No| I
    
    P --> I
    I --> Q{More Cells?}
    Q -->|Yes| F
    Q -->|No| R{Sort Descending?}
    
    R -->|Yes| S[Sort by Value Descending]
    R -->|No| T[Return Local Maxima List]
    S --> T
```

## 5. Grid Operations Workflow

```mermaid
flowchart TD
    A[Grid Operations] --> B{Operation Type?}
    
    B -->|Addition| C[Grid Addition]
    B -->|Coordinate Conversion| D[Coordinate Conversion]
    B -->|Data Access| E[Data Access]
    B -->|Statistics| F[Statistics Computation]
    
    C --> G[Check Grid Dimensions Match]
    G --> H{Match?}
    H -->|No| I[Throw Error]
    H -->|Yes| J[Add Corresponding Cells]
    J --> K[Return New Grid]
    
    D --> L{Conversion Type?}
    L -->|Pixel to Cartesian| M[Calculate Cartesian Coordinates]
    L -->|Cartesian to Pixel| N[Calculate Pixel Indices]
    M --> O[Return 3D Vector]
    N --> P[Return Pixel Coordinates]
    
    E --> Q[Get Raw Data Pointer]
    Q --> R[Return Data Array]
    
    F --> S[Initialize Min/Max]
    S --> T[Scan All Cells]
    T --> U[Update Min/Max Values]
    U --> V[Return Statistics]
```

## 6. Complete Processing Pipeline

```mermaid
flowchart TD
    A[Point Cloud File] --> B[Read Point Cloud Data]
    B --> C[Extract Points and Normals]
    C --> D[STL_STEPCreateGrid3D]
    D --> E[Create 3D Grid]
    E --> F{Apply Filters?}
    
    F -->|Yes| G{Filter Type?}
    F -->|No| H[Final Grid Output]
    
    G -->|Neighbor| I[STL_StepFilterGrid3D]
    G -->|Threshold| J[STL_StepFilterGrid3DByValue]
    G -->|Ratio| K[STL_StepFilterByRatio]
    
    I --> L[Apply Neighbor Filtering]
    J --> M[Apply Threshold Filtering]
    K --> N[Apply Fast Filter]
    
    L --> O[Filtered Grid]
    M --> O
    N --> O
    
    O --> P{Detect Local Maxima?}
    P -->|Yes| Q[Find Local Maxima]
    P -->|No| H
    
    Q --> R[Extract Peak Positions]
    R --> S[Sort by Value]
    S --> T[Final Results]
    
    H --> U[Grid Statistics]
    T --> U
    U --> V[Export Results]
```

## 7. Multi-threading in Grid Creation

```mermaid
flowchart TD
    A[Point Cloud with N Points] --> B[Calculate Thread Count]
    B --> C[Divide Points Among Threads]
    C --> D[Create Thread Pool]
    
    D --> E[Thread 1: Points 0 to N/4]
    D --> F[Thread 2: Points N/4 to N/2]
    D --> G[Thread 3: Points N/2 to 3N/4]
    D --> H[Thread 4: Points 3N/4 to N]
    
    E --> I[Create Local Grid 1]
    F --> J[Create Local Grid 2]
    G --> K[Create Local Grid 3]
    H --> L[Create Local Grid 4]
    
    I --> M[Process Local Points]
    J --> N[Process Local Points]
    K --> O[Process Local Points]
    L --> P[Process Local Points]
    
    M --> Q[Local Grid 1 Complete]
    N --> R[Local Grid 2 Complete]
    O --> S[Local Grid 3 Complete]
    P --> T[Local Grid 4 Complete]
    
    Q --> U[Wait for All Threads]
    R --> U
    S --> U
    T --> U
    
    U --> V[Merge All Grids]
    V --> W[Final Combined Grid]
```

## 8. Beam Visitor Pattern

```mermaid
flowchart TD
    A[Beam Creation] --> B[Initialize Visitor List]
    B --> C[Add STL_Grid3DBeamVisitor]
    C --> D[Add STL_VisitorGrid3DFastFilter]
    D --> E[Add STL_VisitorGrid3DSetValue]
    
    E --> F[Create Woo Traversal Algorithm]
    F --> G[Set Grid and Visitors]
    G --> H[Start Beam Traversal]
    
    H --> I[For Each Cell in Beam Path]
    I --> J[Call Visitor.visit(index, beam)]
    J --> K{Visitor Type?}
    
    K -->|BeamVisitor| L[Add Value to Grid]
    K -->|FastFilter| M[Accumulate Votes]
    K -->|SetValue| N[Set Specific Value]
    
    L --> O[Continue Traversal]
    M --> O
    N --> O
    
    O --> P{More Cells?}
    P -->|Yes| I
    P -->|No| Q[Traversal Complete]
    
    Q --> R[Get Visitor Results]
    R --> S[Process Results]
``` 