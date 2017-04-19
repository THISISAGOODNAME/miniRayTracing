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

> spp同为100的情况下

> 在chrome浏览器上，h5和asm.js两个版本表现十分接近，没有区别

> 在edge浏览器上，asm.js版本渲染完成时h5版本只完成了29%，两者差异巨大

> 结论：如果条件允许，不考虑应用程序容量的情况下，可以把一些重度计算交给asm.js完成

[WebAssembly version](http://aicdg.com/miniRayTracing/smallpt-wasm/smallpt.html?spp=100)

```bash
http://aicdg.com/miniRayTracing/smallpt-wasm/smallpt.html?spp=100
```

> 在chrome浏览器上，WebAssembly版本渲染完成时html5版本只渲染了18%，两者差异异常巨大

> 容量上，asm.js版本最终文件总容量270kb，wasm版本最终总文件容量259kb，wasm版本缩小了4%，虽然不及官方描述的5~7%，但一个个例真心说明不了什么。相比asm.js容量更小是货真价实的。wasm四大主流浏览器都已经开启了支持，但并不能像asm.js一样在浏览器不支持时保持兼容。对于部分项目可以考虑采用

# 其他语言版本

[叶劲峰老师的其他语言fork以及性能比较](http://www.cnblogs.com/miloyip/archive/2010/07/07/languages_brawl_GI.html)

