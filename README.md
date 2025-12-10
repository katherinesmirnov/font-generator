## Font Generator
### About
Takes an input string of text and outputs a B-spline piecewise representation.

### In action:
!["Command line example of 'hello world'"](https://github.com/katherinesmirnov/font-generator/tree/master/photos/cmd_output.png)
!["Command line example of the alphabet"](https://github.com/katherinesmirnov/font-generator/tree/master/photos/alphabet.png)

#### Docker Commands:
```
docker build --tag font-generator
docker run  -it font-generator
docker cp <container ID>:plot.pdf <target file>.pdf

```
