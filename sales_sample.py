#! /usr/bin/env python
#coding=utf-8
import pypyodbc
pypyodbc.win_create_mdb('D:\\salesdb.mdb')
conn = pypyodbc.connect('Driver={Microsoft Access Driver (*.mdb)};DBQ=D:\\salesdb.mdb')
cur = conn.cursor()
cur.execute('''CREATE TABLE saleout (
ID COUNTER PRIMARY KEY,
customer_name VARCHAR(25), 
product_name VARCHAR(30), 
price float, 
volume int,
sell_time datetime);''')

cur.commit()


cur.execute('''INSERT INTO saleout(customer_name,product_name,price,volume,sell_time) 
VALUES(?,?,?,?,?)''',(u'江文','Huawei Ascend mate','5000.5',2,'2012-1-21'))
cur.commit()



cur.execute('''INSERT INTO saleout(customer_name,product_name,price,volume,sell_time) 
VALUES(?,?,?,?,?)''',(u'杨天真','Apple IPhone 5','5500.1',1,'2012-1-21'))
cur.execute('''INSERT INTO saleout(customer_name,product_name,price,volume,sell_time) 
VALUES(?,?,?,?,?)''',(u'郑现实','Huawei Ascend D2','5100.5',1,'2012-1-22'))
cur.execute('''INSERT INTO saleout(customer_name,product_name,price,volume,sell_time) 
VALUES(?,?,?,?,?)''',(u'莫小闵','Huawei Ascend D2','5000.5',1,'2012-1-22'))
cur.execute('''INSERT INTO saleout(customer_name,product_name,price,volume,sell_time) 
VALUES(?,?,?,?,?)''',(u'顾小白','Huawei Ascend mate','5000.5',1,'2012-1-22'))

cur.commit()
cur.execute('''SELECT * FROM saleout WHERE product_name LIKE '%Huawei%';''')
for d in cur.description: print d[0],
print ''
for row in cur.fetchall():
    for field in row: 
        print field,
    print ''
    

conn.close()

pypyodbc.win_compact_mdb('D:\\salesdb.mdb','D:\\salesdb_backup.mdb')
