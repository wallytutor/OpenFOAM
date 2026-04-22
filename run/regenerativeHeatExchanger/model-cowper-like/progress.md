# Progress tracking

```python
import heat_recovery.cowper_like as cl
```

```python
# cl.plot_convergence()
```

## Preliminary

```python
origin = "0.000000e+00"
```

```python
p = cl.plot_temperature("Initialization", origin, loc=3)
```

```python
p = cl.plot_pressure("Initialization", origin)
```

```python
p = cl.plot_flowrate("Initialization", origin, loc=3)
```

```python
p = cl.plot_table("solid", "solidTemperature", "volFieldValue", time=origin)
```

## Relaxation 1

```python
origin = "2.000000e+00"
```

```python
p = cl.plot_temperature("Initialization", origin, loc=3)
```

```python
p = cl.plot_pressure("Initialization", origin)
```

```python
p = cl.plot_flowrate("Initialization", origin, loc=3)
```

```python
p = cl.plot_table("solid", "solidTemperature", "volFieldValue", time=origin)
```

## Discharging

```python
origin = "3.000000e+00"
```

```python
# p = cl.plot_temperature("Initialization", origin, loc=3)
# p = cl.plot_pressure("Initialization", origin)
# p = cl.plot_flowrate("Initialization", origin, loc=3)
# p = cl.plot_table("solid", "solidTemperature", "volFieldValue", time=origin)
```
