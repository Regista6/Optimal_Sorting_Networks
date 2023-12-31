## Optimal_Sorting_Networks

- Model based on [Optimal Sorting Networks (Daniel Bundala and Jakub Zavodny)](https://arxiv.org/pdf/1310.6271.pdf)
- Requires [OR-Tools](https://github.com/google/or-tools) and C++20 (because of `std::format`)

### Example Output
```
Channel Size = 4
Original Depth = 4
Obj_Type = MINIMIZE_DEPTH, Obj_Val = 3
Depth: 1, Index: [(0, 3), (1, 2)]
Depth: 2, Index: [(0, 1), (2, 3)]
Depth: 3, Index: [(1, 2)]
Validation Started
Validation Finished
```
