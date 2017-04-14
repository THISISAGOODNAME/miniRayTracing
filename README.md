# miniRayTracing
编写一个超微型光线追踪渲染器

## c++版

[http://www.kevinbeason.com/smallpt/](http://www.kevinbeason.com/smallpt/)

```bash
# 开启openMP
g++ -O3 -fopenmp smallpt.cpp -o smallpt 
# 不启用openMP
g++ -O3 smallpt.cpp -o smallpt 
time ./smallpt 5000
display image.ppm
```

> C++不要在c11编译，会有int到unsigned short*转换时发生类型窄化的错误

## HTML5版

[Html5 version](http://aicdg.com/miniRayTracing/smallpt.html?spp=100)

```bash
http://aicdg.com/miniRayTracing/smallpt.html?spp=100
```

[asm.js version](http://aicdg.com/miniRayTracing/smallpt-asmjs/smallpt.html?spp=100)

```bash
http://aicdg.com/miniRayTracing/smallpt-asmjs/smallpt.html?spp=100
```

# 其他语言版本

[叶劲峰的其他语言fork以及性能比较](http://www.cnblogs.com/miloyip/archive/2010/07/07/languages_brawl_GI.html)

