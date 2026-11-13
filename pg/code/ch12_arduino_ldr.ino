// ch12_arduino_ldr.ino
//
// 原书 ldr.ino：读 A0/A1 两路 LDR 模拟值，以"值 值\n"格式发串口（每 50 ms 一行）
// 引脚：LDR 接 A0、A1；Serial.begin(9600) 与 Python 端一致

void setup()
{
  // initialize serial comms
  Serial.begin(9600);
}

void loop()
{
  // read A0
  int val1 = analogRead(0);
  // read A1
  int val2 = analogRead(1);
  // print to serial
  Serial.print(val1);
  Serial.print(" ");
  Serial.print(val2);
  Serial.print("\n");
  // wait
  delay(50);
}
