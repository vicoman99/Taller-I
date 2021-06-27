Qué hacer para compilar:

Lo ideal es llamar al proyecto desde CLion, que es donde he 
estado generando el código y donde está probado que funciona.

Sino, desde la consola de comandos, lo que se puede hacer es:

Usar el comando gcc llamando a todos los ".c" al mismo tiempo,
sino da un problema de referencias indefinidas.
El proceso es más sencillo ya que he incluido todos los ficheros
en la misma carpeta, excepto los ficheros de datos con las matrices
del main.

Incluir además "-lm" por el uso de la librería "math.h"

Al tener dos mains (el de los tests y el del principal), yo comento
el que no voy a usar con /*  ...  */ (el documento entero de arriba
a abajo), por ejemplo ahora los tests están comentados y el main 
está descomentado, cuando se quiera ejecutar tests, comentar el main y
descomentar los tests.