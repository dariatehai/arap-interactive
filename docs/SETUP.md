<!-- Start commands -->

## Build
```
cd /globalpath/to/arap-interactive
mkdir -p build && cd build
cmake .. && cmake --build . && ./tests/arap_tests
```


## For debugging
```
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
cmake --build build --config Debug --target arap_app 
.\build\src\Debug\arap_app.exe 
```


## Links to models
[Simple PLYs](hhttps://people.sc.fsu.edu/~jburkardt/data/ply/ply.html).


## Github commands
If you need to clone:
```
git clone URL
```

Pulling from the main branch
```
git pull origin main
```

Add the files before commit (the ones you changed)
```
git add .
```

```
git commit -m "your message"
```

To push exclusively to your branch
```
git switch -c your_branch_name
git push origin your_branch_name
```